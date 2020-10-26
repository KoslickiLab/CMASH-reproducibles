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
import screed
import bisect

#####################
# local variables for convenience
num_threads = 6
prime = 9999999999971
ksize = 60
input_range = "10-60-5"
max_h = 1000
rev_comp = True
fast_mode = False
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

# functions copied from CMash scripts
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

def is_prime(number):
	if not isinstance(number, int):
		raise Exception("Input number is not an integer")
	if number < 2:
		return False
	if number == 2:
		return True
	if number % 2 == 0:
		return False
	for _ in range(3, int(number ** 0.5) + 1, 2):
		if number % _ == 0:
			return False
	return True

def get_prime_lt_x(input_value):
	"""
	Backward-find a prime smaller than the input value
	:param input_value: any number >=2
	:return: a primer number <= input
	"""
	if input_value < 1:
		raise Exception("Input value <1!")
	if input_value == 1:
		return 1
	i = int(input_value)
	if i % 2 == 0:
		i -= 1
	while i > 0:
		if is_prime(i):
			return i
		i -= 2

def kmers(seq, ksize):
	for i in range(len(seq) - ksize + 1):
		yield seq[i:i+ksize]
	
# wrapper for trunc_JI from JI_CE object
def unwrap_jaccard_vector(arg):
		return arg[0].est_jaccard(arg[1])



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



# create kmc database for all files for all k values
def create_kmc_db(input_file, k):
	out_name = f"{os.path.join(temp_dir, os.path.basename(input_file))}_k_"
	out_file_name = out_name + str(k)
	res = subprocess.run(
		f"kmc -k{k} -fm -r -t{num_threads} -cs1 -ci0 -j{out_file_name}.log {input_file} {out_file_name} ."
		, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
	if res.returncode != 0:
		raise Exception(f"The command {res.args} failed to run and returned the returncode={res.returncode}")

if not fast_mode:
	print("Creating hmc database for JI test")
	if not os.path.exists(temp_dir):
		os.mkdir(temp_dir)
	for file in set(query_list).union(set(ref_list)):
		for k in k_sizes:
			print("Running file " + file)
			create_kmc_db(file, k)



# get ground truth JI by kmc
def kmc_intersection_count_for_ji(kmc_input_file1: str, kmc_input_file2: str, input_k, keep_intermediate=True):
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
	res = subprocess.run(f"kmc_dump temp_intersect_min temp_dump_intersect; cat temp_dump_intersect | wc -l | cut -f 1", shell=True,
						 capture_output=True)
	ji_numerator = int(res.stdout)
	print(ji_numerator)
	res = subprocess.run(f"kmc_dump temp_union_min temp_dump_union; cat temp_dump_union | wc -l | cut -f 1",
						 shell=True, capture_output=True)
	ji_denominator = int(res.stdout)
	print(ji_denominator)
	# rm intermediate files
	os.remove("temp_intersect_min.kmc_pre")
	os.remove("temp_intersect_min.kmc_suf")
	os.remove("temp_union_min.kmc_pre")
	os.remove("temp_union_min.kmc_suf")
	if keep_intermediate:
		os.chdir("../")
		subprocess.run(f"mv {temp_name} ./{temp_dir}", shell=True )
	else:
		os.remove("temp_dump_union")
		os.remove("temp_dump_intersect")
		os.chdir("../")
		try:
			os.rmdir(temp_name)
		except OSError:
			print("Processing unsucceful for " + temp_name)
			raise SystemExit
	# return JI
	return ji_numerator * 1.0 / ji_denominator

# wrapper for GT JI by kmc
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
### CE object similar to CMash
class JI_CountEstimator(object):
	def __init__(self, n=None, max_prime=9999999999971, ksize=None, input_file_name=None, rev_comp=True):
		if n is None:
			raise Exception
		if ksize is None:
			raise Exception

		# read parameters
		self.ksize = ksize
		self.p = get_prime_lt_x(max_prime)
		self.input_file_name = input_file_name
		self.rev_comp = rev_comp

		# initialize sketch of size n with count
		self._mins = [self.p] * n
		self._kmers = [''] * n

		# read sequence by "screed" and load them into the sketch (replace parse_file function in L105)
		if self.input_file_name:
			for record in screed.open(self.input_file_name):
				self.add_sequence(record.sequence)

	def add_sequence(self, seq):
		seq = seq.upper() #transfer to upper case (just in case)
		seq_split_onlyACTG = re.compile('[^ACTG]').split(seq) # get rid of non-ATCG letter
		if len(seq_split_onlyACTG) == 1:
			if len(seq) >= self.ksize:
				for kmer in kmers(seq, self.ksize):
					self.add(kmer)
		else:
			for sub_seq in seq_split_onlyACTG:
				if len(sub_seq)>=self.ksize:        #in case of small chunk
					self.add_sequence(sub_seq)

	def add(self, kmer):
		_mins = self._mins
		_kmers = self._kmers
		# use rev_comp if needed
		if self.rev_comp:
			cano_kmer = min(kmer, khmer.reverse_complement(kmer))
			h = khmer.hash_no_rc_murmur3(cano_kmer)
		else:
			h = khmer.hash_no_rc_murmur3(kmer)
		# reminder of max_prime we use
		h = h % self.p
		# early stop if n sketches are found
		if h >= _mins[-1]:
			return
		# insert kmer into the sketch
		i = bisect.bisect_left(_mins, h)  # find index to insert h
		if _mins[i] == h: #already in sketch
			return
		else:
			#h not in sketch, insert
			_mins.insert(i, h)
			_kmers.insert(i, kmer)
			_mins.pop()
			_kmers.pop()
			return

	def est_jaccard(self, other):
		if self.ksize != other.ksize:
			raise Exception("different k-mer sizes - cannot compare")
		if self.p != other.p:
			raise Exception("different primes - cannot compare")
		i = 0
		j = 0
		overlap = 0
		processed = 0
		mins1 = self._mins
		mins2 = other._mins
		while processed < min(len(mins1), len(mins2)):  # loop stop when min(|A|, |B|) sample size reached
			if mins1[i] < mins2[j]:
				i += 1
				processed += 1
			elif mins1[i] > mins2[j]:
				j += 1
				processed += 1
			elif mins1[i] == mins2[j]:
				# a match from A overlap B
				overlap += 1
				i += 1
				j += 1
				processed += 1
		print(overlap)
		print(processed)
		est_ji = overlap * 1.0 / processed
		return est_ji

	def ji_vector(self, other_list, threads=num_threads):
		pool = Pool(processes=threads)
		Y = np.array(pool.map(unwrap_jaccard_vector, zip([self] * len(other_list), other_list)))
		pool.terminate()
		return Y

	def brute_force_truncation(self, new_ksize):
		if not isinstance(new_ksize, int):
			raise Exception("Input number is not an integer")
		if new_ksize > self.ksize:
			raise Exception("New size must be smaller than %d." %self.ksize)
		elif new_ksize == self.ksize:
			return
		elif new_ksize < self.ksize:
			# data to be updated after the truncation:
			self.ksize = new_ksize
			while self._mins[-1] == self.p: # rm unused cells, otherwise empty cell (though very rare) has hash value 0
				self._mins.pop()
				self._kmers.pop()
			new_kmers = list(set([x[0:new_ksize] for x in self._kmers]))
			sketch_size = len(new_kmers)
			self._mins = [self.p] * sketch_size
			self._kmers = [''] * sketch_size
			# update
			for i in range(sketch_size):
				self.add(new_kmers[i])
			while self._mins[-1] == self.p: # rm unused cells
				self._mins.pop()
				self._kmers.pop()
			return


### get JI_CE object lists and array functions
def make_minhash(genome, max_h, prime, ksize, rev_comp):
	MHS = JI_CountEstimator(n=max_h, max_prime=prime, ksize=ksize, input_file_name=genome, rev_comp=rev_comp)
	return MHS

def make_minhash_star(arg):
	return make_minhash(*arg)

def get_ce_lists(input_file, input_k, reverse_comp):
	para = Pool(processes=num_threads)
	genome_sketches = para.map(make_minhash_star, zip(input_file, repeat(max_h), repeat(prime), repeat(input_k), repeat(reverse_comp)))
	para.close()
	return genome_sketches

def ji_matrix(list1, list2, out_file_name):
	lst = []
	row_name = []
	col_name = [os.path.basename(list2[i].input_file_name) for i in range(len(list2))]
	for i in range(len(list1)):
		name = os.path.basename(list1[i].input_file_name)
		row_name.append(name)
		Y = list1[i].ji_vector(list2, num_threads)
		lst.append(Y)
	df = pd.DataFrame(lst)
	df.columns = col_name
	df.index = row_name
	df.to_csv(out_file_name, index=True, encoding='utf-8')
	return df

### get trunc_JI
sketch1 = get_ce_lists(query_list, ksize, rev_comp)
sketch2 = get_ce_lists(ref_list, ksize, rev_comp)

def bf_truncation(ce_list, new_ksize):
	for ce in ce_list:
		ce.brute_force_truncation(new_ksize)

for i in k_sizes:
	bf_truncation(sketch1, i)
	bf_truncation(sketch2, i)
	out_name="trunc_JI_k"+str(i)+".csv"
	out1 = ji_matrix(sketch1, sketch2, out_name)

### get est_JI
for i in k_sizes:
	sketch1 = get_ce_lists(query_list, i, rev_comp)
	sketch2 = get_ce_lists(ref_list, i, rev_comp)
	out_name = "est_JI_k"+str(i)+".csv"
	out2 = ji_matrix(sketch1, sketch2, out_name)



# get the bias factor: work on full kmer set (not the random sample) by kmc
### to calculate the bias factor, we need 4 things:
### 1. | A overlap B |
### 2. | A union B |
### 3. E(RE on A/capB)
### 4. E(RE on A/cupB)
### 1,2 are already there from the GT_JI calculation by kmc (keep intermediate!!!)
### 3,4 can be obtained by grepping kmers of k_small from kmers of k_large using prefix: grep ^
def get_bias_factor(file1, file2, maxk, truncated_k, keep_intermediate = False):
	"""
	:param file1: 1 query file
	:param file2: 1 ref file
	:param maxk: full kmer length (usually 60)
	:param truncated_k: truncated kmer length (e.g. 60 -> 55)
	:return: bias factor from maxk -> truncated_k between the 2 input files
	"""
	file1 = os.path.basename(file1)
	file2 = os.path.basename(file2)
	folder_k_small = f"{temp_dir}/kmc_output_{file1}_k_{truncated_k}_against_{file2}_k_{truncated_k}"
	folder_k_large = f"{temp_dir}/kmc_output_{file1}_k_{maxk}_against_{file2}_k_{maxk}"
	
	### get cardinality
	res = subprocess.run(f"cat {folder_k_small}/temp_dump_intersect | wc -l | cut -f 1", shell=True, capture_output=True)
	cardi_ksmall_cap = int(res.stdout)
	res = subprocess.run(f"cat {folder_k_small}/temp_dump_union | wc -l | cut -f 1", shell=True, capture_output=True)
	cardi_ksmall_cup = int(res.stdout)
	
	### get right extension count from k_large union set
	### RE input of ksmall cap
	res = subprocess.run(f"cut -f 1 {folder_k_small}/temp_dump_intersect | awk \'{{ print \"^\"$1 }}\' >> temp_re_input_intersect_{file1}_{file2}_k{truncated_k}.txt", shell=True)
	### RE input of ksmall cup
	res = subprocess.run(f"cut -f 1 {folder_k_small}/temp_dump_union | awk \'{{ print \"^\"$1 }}\' >> temp_re_input_union_{file1}_{file2}_k{truncated_k}.txt",shell=True)
	
	### RE count of ksmall cap
	res = subprocess.run(f"grep -f temp_re_input_intersect_{file1}_{file2}_k{truncated_k}.txt {folder_k_large}/temp_dump_union | wc -l | cut -f 1", shell=True, capture_output=True)
	RE_cap = int(res.stdout)
	### RE count of ksmall cup
	res = subprocess.run(f"grep -f temp_re_input_union_{file1}_{file2}_k{truncated_k}.txt {folder_k_large}/temp_dump_union | wc -l | cut -f 1", shell=True, capture_output=True)
	RE_cup = int(res.stdout)
	
	### bias factor
	bias_calculation = (RE_cap * 1.0 / cardi_ksmall_cap) / (RE_cup * 1.0 / cardi_ksmall_cup)
	
	### clean files
	if not keep_intermediate:
		os.remove(f"temp_re_input_intersect_{file1}_{file2}_k{truncated_k}.txt")
		os.remove(f"temp_re_input_union_{file1}_{file2}_k{truncated_k}.txt")
	
	print(bias_calculation)
	return bias_calculation

def unwrap_bias_factor_vector(arg):
	return get_bias_factor(arg[0], arg[1], arg[2], arg[3])

def get_bias_factor_vector(input_file, other_list, maxk, truncated_k):
	pool = Pool(processes=num_threads)
	Y = np.array(pool.map(unwrap_bias_factor_vector, zip([input_file] * len(other_list), other_list, [maxk] * len(other_list), [truncated_k] * len(other_list))))
	pool.terminate()
	return Y

def bias_factor_matrix(list1, list2, input_k, out_file_name):
	lst = []
	row_name = []
	col_name = [os.path.basename(list2[i]) for i in range(len(list2))]
	for i in range(len(list1)):
		name = os.path.basename(list1[i])
		row_name.append(name)
		Y = get_bias_factor_vector(list1[i], list2, ksize, input_k)     # ksize is the input maxk variable
		lst.append(Y)
	df = pd.DataFrame(lst)
	df.columns = col_name
	df.index = row_name
	df.to_csv(out_file_name, index=True, encoding='utf-8')
	return df

for k in k_sizes:
	out_name = "bias_factor_of_k"+str(k)+"_against_k"+str(ksize)+".csv"
	bias_factor_matrix(query_list, ref_list, k, out_name)




##########################################################################
# Tests
### test ground truth JI
def f1_test_gt_ji():
	subprocess.run(f"echo TTAATTAA > test1.fq", shell=True)
	subprocess.run(f"echo AAATTTAA > test2.fq", shell=True)
	create_kmc_db("test1.fq", 4)
	create_kmc_db("test2.fq", 4)
	test_ji = kmc_intersection_count_for_ji("test1.fq", "test2.fq", 4)
	os.remove("test1.fq")
	os.remove("test2.fq")
	print("This is for canonical kmers: |A|=3, |B|=4, overlap=2")
	print(test_ji)
	assert test_ji == 2 / 5.0

### test est JI
def f2_test_est_ji():
	E1 = JI_CountEstimator(n=0, ksize=21)
	E2 = JI_CountEstimator(n=0, ksize=21)
	E1._mins = [1, 2, 3, 4, 5]
	E2._mins = [1, 2, 3, 4]
	assert E1.est_jaccard(E2) == 4/4.0
	assert E2.est_jaccard(E1) == 4/4.0
	E2._mins = [1, 3, 4, 5, 6]
	assert E1.est_jaccard(E2) == 4/5.0
	
### test truncation
def f3_test_truncation():
	E1 = JI_CountEstimator(n=5, ksize=4)
	E2 = JI_CountEstimator(n=5, ksize=4)
	E1.add_sequence("TTAATTAA")
	print(E1._kmers)
	E1.brute_force_truncation(3)
	print(E1._kmers)
	assert E1._kmers == ['AAT', 'TTA']
	