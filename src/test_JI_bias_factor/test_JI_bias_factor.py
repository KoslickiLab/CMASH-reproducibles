# read parameters
# This is for unweighted JI (not CI), to test the bias factor
# 2 required: -q <query_file> -r <ref_file>
import multiprocessing
from multiprocessing import Pool
import argparse
import os
import sys
import khmer
from itertools import *
import numpy as np
import pandas as pd
import re
from argparse import ArgumentTypeError
from CMash import MinHash as MH
from CMash.Make import MakeTSTNew
import marisa_trie as mt
import subprocess

#####################
# local variables for convenience
num_threads = 6
prime = 9999999999971
ksize = 60
k_sizes = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
max_h = 1000
rev_comp = True
query_file = "/Users/shaopeng/Desktop/WTST_test_run/query_path.txt"
ref_file = "/Users/shaopeng/Desktop/WTST_test_run/ref_path.txt"
#####################




parser = argparse.ArgumentParser(description="This script creates training/reference sketches for each FASTA/Q file"
									" listed in the input file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-q','--query', help="Path to the query file.")
parser.add_argument('-r', '--ref' , help="Path to the ref file.")
parser.add_argument('-f', '--fast' , type=str, help="Skipping steps to build GT_TST and CE.h5 file", default="False")
parser.add_argument('-c', '--rev_comp' , type=str, help="Whether to keep the reverse complementary", default="True")
parser.add_argument('-n', '--num_hashes', type=int, help="Number of hashes to use.", default=1000)
parser.add_argument('-k', '--k_size', type=int, help="k-mer size", default=60)
parser.add_argument('-g', '--k_range', type=str, help="k-mer range", default="10-60-5")
parser.add_argument('-t', '--threads', type=int, help="Number of threads to use", default=min(64, int(multiprocessing.cpu_count()/2)))
parser.add_argument('-p', '--prime', help='Prime (for modding hashes)', default=9999999999971)

args = parser.parse_args()
num_threads = args.threads
prime = args.prime  # taking hashes mod this prime
ksize = args.k_size # max k used in the analysis
max_h = args.num_hashes # number of hashes to use
input_range = args.k_range # k range
query_file = args.query
ref_file = args.ref
rev_comp = args.rev_comp
fast_mode = args.fast
if rev_comp == "True":
	rev_comp = True
elif rev_comp == "False":
	rev_comp = False
if rev_comp:
	print("Using true for rev_comp!!!!!!")

if fast_mode == "False":
	fast_mode = False
elif fast_mode == "True":
	fast_mode = True
if fast_mode:
	print("Skipping all kmc building process and assume they are already built!!!!!!")

query_file_names = os.path.abspath(query_file)
if not os.path.exists(query_file_names):
	raise Exception("Input file %s does not exist." % query_file_names)
ref_file_names = os.path.abspath(ref_file)
if not os.path.exists(ref_file_names):
	raise Exception("Input file %s does not exist." % query_file_names)

def parsenumlist(k_sizes_str: str):
	"""
	Parses a string like 10-21-1 and turn it into a list like [10, 11, 12,...,21]
	:param k_sizes_str: the <start>-<end>-<increment> string
	:type k_sizes_str: str
	:return: list of k-mer sizes
	:rtype: list
	"""
	m = re.match(r'(\d+)(?:-(\d+))?(?:-(\d+))?$', k_sizes_str)
	# ^ (or use .split('-'). anyway you like.)
	if not m:
		raise ArgumentTypeError(
			"'" + k_sizes_str + "' is not a range of number. Expected forms like '1-5' or '2' or '10-15-2'.")
	start = int(m.group(1))
	end = int(m.group(2))
	if m.group(3):
		increment = int(m.group(3))
	else:
		increment = 1
	return list(range(start, end + 1, increment))
k_sizes=parsenumlist(input_range)

# functions copied from CMash scripts
def make_minhash(genome, max_h, prime, ksize):
	MHS = MH.CountEstimator(n=max_h, max_prime=prime, ksize=ksize, save_kmers='y', input_file_name=genome, rev_comp=False)  # the query automatically takes care of rev_comp's for me
	# Just use HLL to estimate the number of kmers, no need to get exact count
	hll = khmer.HLLCounter(0.01, ksize)
	hll.consume_seqfile(genome)
	MHS._true_num_kmers = hll.estimate_cardinality()
	MHS.input_file_name = genome
	return MHS

def make_minhash_star(arg):
	return make_minhash(*arg)


# pipe start
def check_files(input_file):
	print("Reading paths from %s to a list." %input_file)
	out_list = list()
	temp = open(input_file, 'r')
	for line in temp.readlines():
		line = line.strip()
		if not os.path.exists(line):
			raise Exception("Input file %s does not exist." % line)
		out_list.append(os.path.abspath(line))
	temp.close()
	out_list = sorted(out_list, key=os.path.basename)
	return(out_list)
query_list = check_files(query_file)
ref_list = check_files(ref_file)
temp_dir = "JI_test_kmc_database"

# while the JI is symmatric, this time I only check query against all refs
def create_kmc_db(input_file):
	out_name = f"{os.path.join(temp_dir, os.path.basename(input_file))}_k_"
	for k in k_sizes:
		out_file_name = out_name + str(k)
		res = subprocess.run(
			f"kmc -k{k} -fm -r -t{num_threads} -cs1 -ci0 -j{out_file_name}.log {input_file} {out_file_name} ."
			, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
		if res.returncode != 0:
			raise Exception(f"The command {res.args} failed to run and returned the returncode={res.returncode}")

# create kmc kmer files for all input data
if not fast_mode:
	print("Creating hmc database for JI test")
	if not os.path.exists(temp_dir):
		os.mkdir(temp_dir)
	for file in set(query_list).union(set(ref_list)):
		print("Running file " + file)
		create_kmc_db(file)



# get ground truth JI
def kmc_intersection_count_for_ji(kmc_input_file1: str, kmc_input_file2: str, input_k, keep_intermediate=False):
	# use intermediate folder to avoid file name conflicts
	kmc_input_file1 = os.path.basename(kmc_input_file1)
	kmc_input_file2 = os.path.basename(kmc_input_file2)
	temp_name = f"kmc_output_{kmc_input_file1}_k_{input_k}_against_{kmc_input_file2}_k_{input_k}"
	os.mkdir(temp_name)
	os.chdir(temp_name)
	input1 = "../" + f"{os.path.join(temp_dir, kmc_input_file1)}" + "_k_" + str(input_k)
	input2 = "../" + f"{os.path.join(temp_dir, kmc_input_file2)}" + "_k_" + str(input_k)
	# numerator: intersect -ocmin
	# denominator: union -ocmin
	res = subprocess.run(f"kmc_tools simple {input1} -ci1 {input2} -ci1 intersect temp_intersect_min -ocmin union temp_union_min -ocmin",
						shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	if res.returncode != 0:
		raise Exception(f"The command {res.args} failed to run and returned the returncode={res.returncode}")
	# get results
	res = subprocess.run(f"kmc_dump temp_intersect_min temp_dump_intersect; cat temp_dump_intersect | wc -l | cut -f 2", shell=True,
						 capture_output=True)
	ji_numerator = int(res.stdout)
	print(ji_numerator)
	res = subprocess.run(f"kmc_dump temp_union_min temp_dump_union; cat temp_dump_union | wc -l | cut -f 2",
						 shell=True, capture_output=True)
	ji_denominator = int(res.stdout)
	print(ji_denominator)
	# rm intermediate files
	os.remove("temp_intersect_min.kmc_pre")
	os.remove("temp_intersect_min.kmc_suf")
	os.remove("temp_union_min.kmc_pre")
	os.remove("temp_union_min.kmc_suf")
	if not keep_intermediate:
		os.remove("temp_dump_union")
		os.remove("temp_dump_intersect")
		os.chdir("../")
		try:
			os.rmdir(temp_name)
		except OSError:
			print("Processing unsucceful for " + temp_name)
			raise SystemExit
	# return WJI
	return ji_numerator * 1.0 / ji_denominator

def unwrap_ji_vector(arg):
	return kmc_intersection_count_for_ji(arg[0], arg[1], arg[2])

def get_ji_vector(input_file, other_list, input_k, threads=num_threads):
	pool = Pool(processes=threads)
	Y = np.array(pool.map(unwrap_ji_vector, zip([input_file] * len(other_list), other_list, [input_k] * len(other_list))))
	pool.terminate()
	return Y

### generate GT JI matrix
def ji_matrix(list1, list2, input_k, out_file_name):
	lst = []
	row_name = []
	col_name = [os.path.basename(list2[i]) for i in range(len(list2))]
	for i in range(len(list1)):
		name = os.path.basename(list1[i])
		row_name.append(name)
		Y = get_ji_vector(list1[i], list2, input_k)
		lst.append(Y)
	df = pd.DataFrame(lst)
	df.columns = col_name
	df.index = row_name
	df.to_csv(out_file_name, index=True, encoding='utf-8')
	return df

for k in k_sizes:
	out_name = "GroundTruth_JI_k"+str(k)+".csv"
	ji_matrix(query_list, ref_list, k, out_name)



# get estimated JI and truncated JI
### use 1000 minhash values (max_h = 1000) to select kmers and find the overlap in python
### iteratively truncate for subsets
# still use max(|A|, |B|) during truncation?




# get the bias factor: work on full kmer set (not the random sample)
### BF: create kmer set of maxk -> for each smaller k, build a dict with counts (store those dict to save MEM if needed)
###                             ->
### kmc for prefix match search?
### use trie?
		

