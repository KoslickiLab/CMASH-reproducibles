"""
Read 2 files of query/ref, calculate GT_WJI
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
parser.add_argument('-c', '--rev_comp' , type=str, help="Whether to keep the reverse complementary", default="True")
parser.add_argument('-k', '--k_range', type=str, help="k-mer size", default="10-60-5")
parser.add_argument('-t', '--threads', type=int, help="Number of threads to use", default=min(64, int(multiprocessing.cpu_count()/2)))

args = parser.parse_args()
num_threads = args.threads
input_range = args.k_range # max k used in the analysis
query_file = args.query
ref_file = args.ref
rev_comp = args.rev_comp
if rev_comp == "True":
	rev_comp = True
elif rev_comp == "False":
	rev_comp = False
	
query_file_names = os.path.abspath(query_file)
if not os.path.exists(query_file_names):
	raise Exception("Input file %s does not exist." % query_file_names)
ref_file_names = os.path.abspath(ref_file)
if not os.path.exists(ref_file_names):
	raise Exception("Input file %s does not exist." % query_file_names)

### manual test parameters:
# input_range="10-60-5"
# query_file = "/Users/shaopeng/Desktop/WTST_test_run/query_path.txt"
# ref_file = "/Users/shaopeng/Desktop/WTST_test_run/ref_path.txt"
# rev_comp = False
# num_threads = 6
#############


query_list = fc.check_files(query_file)
ref_list = fc.check_files(ref_file)
k_range = fc.parsenumlist(input_range)

# test ground truth
# 1. build list of GT objects
def make_gt(seq_file=None, k_range=None,  rev_comp=None, temp_dir=None, threads=6):
	# a wrapper, may add additional tags
	GTS = fc.GroundTruth_WJI(k_range=k_range, seq_file=seq_file, rev_comp=rev_comp, temp_dir=temp_dir, threads=threads)
	return GTS

def make_gt_star(arg):
	return make_gt(*arg)

def get_gt_list(input_file, input_range, temp_dir, rev_comp, thread_num):
	para = Pool(processes=thread_num)
	GT_obj = para.map(make_gt_star, zip(input_file, repeat(input_range), repeat(rev_comp), repeat(temp_dir), repeat(thread_num)))
	para.close()
	return GT_obj


# this is a wired situation: (hypothesis) in parallel creating GT objects, several if check will be met simoutaneously
# however when call "mkdir", only first will success and others return errors "file exists"
# if I manually create this folder, then there is no error.

if not os.path.exists("kmc_database"):
	os.mkdir("kmc_database")
GT1 = get_gt_list(query_list, input_range, "kmc_database", rev_comp, num_threads)
GT2 = get_gt_list(ref_list, input_range, "kmc_database", rev_comp, num_threads)

# 2. generate all kmer database for k range
def generate_kmc_db(gt_list):
	for gt in gt_list:
		gt._compute_all_training_kmers()

generate_kmc_db(GT1)
generate_kmc_db(GT2)

# 3. generate GT WJI for 1 vs array
def gt_wji_matrix(list1, list2, k_size, out_file_name):
	lst = []
	row_name = []
	col_name = [os.path.basename(list2[i].seq_file) for i in range(len(list2))]
	for i in range(len(list1)):
		name = os.path.basename(list1[i].seq_file)
		row_name.append(name)
		Y = list1[i].gt_wji_vector(list2, k_size)   #thread is in gt_wji function
		lst.append(Y)
	df = pd.DataFrame(lst)
	df.columns = col_name
	df.index = row_name
	df.to_csv(out_file_name, index=True, encoding='utf-8')
	return df

for k_size in k_range:
	out_name = "GT_WJI_k" + str(k_size) + ".csv"
	gt_wji_matrix(GT1, GT2, k_size, out_name)
	
