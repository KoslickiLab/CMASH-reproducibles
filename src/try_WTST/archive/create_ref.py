"""
Read 2 files of query/ref, build WTST for each of them
"""


### read parameter
# 2 required: -q <query_file> -r <ref_file>
import multiprocessing
import argparse
import os
import py_src as fc
import khmer

parser = argparse.ArgumentParser(description="This script creates training/reference sketches for each FASTA/Q file"
									" listed in the input file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-q','--query', help="Path to the query file.")
parser.add_argument('-r', '--ref' , help="Path to the ref file.")
parser.add_argument('-n', '--num_hashes', type=int, help="Number of hashes to use.", default=2000)
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



### start!
#############################################
### Manually input parameters (rm this section after test)
max_h = 500
num_threads = 8
ksize = 40
prime = 9999999999971
query_file = "/Users/shaopeng/Desktop/WTST_test_run/query_path.txt"
ref_file = "/Users/shaopeng/Desktop/WTST_test_run/ref_path.txt"
#############################################

# read each single files into a list
query_list = fc.check_files(query_file)
ref_list = fc.check_files(ref_file)

# go to workdir
if not os.path.exists('test_create_ref'):
	os.makedirs('test_create_ref')
os.chdir("./test_create_ref")

# single file test
f1 = fc.CountEstimator(n=max_h, max_prime=prime, ksize=ksize, input_file_name=ref_list[0])
f2 = fc.CountEstimator(n=max_h, max_prime=prime, ksize=ksize, input_file_name=query_list[0])
f1.est_weighted_jaccard(f2)


