"""
Read 2 files of query/ref, calculate est_WJI and truncated_WJI
"""


### read parameter
# 2 required: -q <query_file> -r <ref_file>
import multiprocessing
from multiprocessing import Pool
import argparse
import os
import sys
import py_src as fc
import khmer
from itertools import *
import numpy as np
import pandas as pd


parser = argparse.ArgumentParser(description="This script creates training/reference sketches for each FASTA/Q file"
									" listed in the input file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-q','--query', help="Path to the query file.")
parser.add_argument('-r', '--ref' , help="Path to the ref file.")
parser.add_argument('-c', '--rev_comp' , help="Whether to keep the reverse complementary", default=True)
parser.add_argument('-n', '--num_hashes', type=int, help="Number of hashes to use.", default=1000)
parser.add_argument('-k', '--k_size', type=int, help="k-mer size", default=60)
parser.add_argument('-t', '--threads', type=int, help="Number of threads to use", default=min(8, int(multiprocessing.cpu_count()/2)))
parser.add_argument('-p', '--prime', help='Prime (for modding hashes)', default=9999999999971)

args = parser.parse_args()
num_threads = args.threads
prime = args.prime  # taking hashes mod this prime
ksize = args.k_size # max k used in the analysis
max_h = args.num_hashes # number of hashes to use
query_file = args.query
ref_file = args.ref
rev_comp = args.rev_comp
query_file_names = os.path.abspath(query_file)
if not os.path.exists(query_file_names):
	raise Exception("Input file %s does not exist." % query_file_names)
ref_file_names = os.path.abspath(ref_file)
if not os.path.exists(ref_file_names):
	raise Exception("Input file %s does not exist." % query_file_names)

############################################
### test parameter input and modules (rm this section after test)
fc.test_run() #test import
print(os.getcwd())
print(query_file)
print("The ref file is " + ref_file)
#############################################


#############################################
### Manually input parameters (rm this section after test)
# max_h = 500
# num_threads = 8
# ksize = 40
# prime = 9999999999971
# threads = 6
# query_file = "/Users/shaopeng/Desktop/WTST_test_run/query_path.txt"
# ref_file = "/Users/shaopeng/Desktop/WTST_test_run/ref_path.txt"
#
# # read each single files into a list
# query_list = fc.check_files(query_file)
# ref_list = fc.check_files(ref_file)
#
# # single file test
# f1 = fc.CountEstimator(n=max_h, max_prime=prime, ksize=ksize, input_file_name=ref_list[0])
# f2 = fc.CountEstimator(n=max_h, max_prime=prime, ksize=ksize, input_file_name=query_list[0])
# f1.est_weighted_jaccard(f2) #0.791
#############################################

# pipe start
query_list = fc.check_files(query_file)
ref_list = fc.check_files(ref_file)

# parellel: build CE objects for both query and ref files
def make_minhash(genome, max_h, prime, ksize, reverse_comp=True):
	# a wrapper, may add additional tags
	MHS = fc.CountEstimator(n=max_h, max_prime=prime, ksize=ksize, input_file_name=genome, rev_comp=reverse_comp)
	return MHS

# add a tup wrapper to correctly recognize input parameters:
# see: https://stackoverflow.com/questions/25075809/multiprocessing-pool-map-function-issue
def make_minhash_star(arg):
	return make_minhash(*arg)

# use pool.imap to make an iteratable object to export
def get_ce_lists(input_file, thread_num, input_k, reverse_comp):
	para = Pool(processes=thread_num)
	genome_sketches = para.map(make_minhash_star, zip(input_file, repeat(max_h), repeat(prime), repeat(input_k), repeat(reverse_comp)))
	para.close()
	return genome_sketches

sketch1 = get_ce_lists(query_list, num_threads, ksize, rev_comp)
sketch2 = get_ce_lists(ref_list, num_threads, ksize, rev_comp)

# store those data into a hdf5 file


#





#### below should be put into another py scipt as running them
# maxk=60
# range=10-60-5

# WJI vector to array
def wji_array(list1, list2, thread_num, out_file_name):
	lst = []
	row_name = []
	col_name = [os.path.basename(list2[i].input_file_name) for i in range(len(list2))]
	for i in range(len(list1)):
		name = os.path.basename(list1[i].input_file_name)
		row_name.append(name)
		Y = list1[i].wji_vector(list2, thread_num)
		lst.append(Y)
	df = pd.DataFrame(lst)
	df.columns = col_name
	df.index = row_name
	df.to_csv(out_file_name, index=True, encoding='utf-8')
	return df

# usually use maxk=60
out_name= "est_WJI_k"+str(ksize)+".csv"
out1 = wji_array(sketch1, sketch2, num_threads, out_name)

# truncation (no need to parallel, because the sketch is small)
def bf_truncation(ce_list, new_ksize):
	for ce in ce_list:
		ce.brute_force_truncation(new_ksize)

truncation_seq = list(range(20, 60, 5))
truncation_seq.reverse()
for i in truncation_seq:
	bf_truncation(sketch1, i)
	bf_truncation(sketch2, i)
	out_name="trunc_WJI_k"+str(i)+".csv"
	out1 = wji_array(sketch1, sketch2, num_threads, out_name)


# standard est_wji data
for i in list(range(20, 60, 5)):
	sketch1 = get_ce_lists(query_list, num_threads, i)
	sketch2 = get_ce_lists(ref_list, num_threads, i)
	out_name = "est_WJI_k" + str(i) + ".csv"
	out1 = wji_array(sketch1, sketch2, num_threads, out_name)









