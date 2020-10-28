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
import subprocess
import screed
import bisect
import time
# sys.path.append('/Users/shaopeng/Desktop/WTST_test_run/test_JI_bias_factor')

#####################
# local variables for convenience
# num_threads = 6
# prime = 9999999999971
# ksize = 60
# k_sizes = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
# input_range = "10-60-5"
# max_h = 1000
# rev_comp = True
# fast_mode = False
# query_file = "/Users/shaopeng/Desktop/WTST_test_run/query_path.txt"
# ref_file = "/Users/shaopeng/Desktop/WTST_test_run/ref_path.txt"
#####################


### function defi

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

# generate kmers of ksize for a sequence
def kmers(seq, ksize):
	for i in range(len(seq) - ksize + 1):
		yield seq[i:i + ksize]

# transform file into a list
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



### GT JI
# create kmc database
def create_kmc_db(input_file, k, thread=6):
	temp_dir = "JI_test_kmc_database"
	if not os.path.exists(temp_dir):
		os.mkdir(temp_dir)
	out_name = f"{os.path.join(temp_dir, os.path.basename(input_file))}_k_"
	out_file_name = out_name + str(k)
	res = subprocess.run(
		f"kmc -k{k} -fm -r -t{thread} -cs1 -ci0 -j{out_file_name}.log {input_file} {out_file_name} ."
		, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
	if res.returncode != 0:
		raise Exception(f"The command {res.args} failed to run and returned the returncode={res.returncode}")
# get ground truth JI by kmc
def kmc_intersection_count_for_ji(kmc_input_file1: str, kmc_input_file2: str, input_k, keep_intermediate=False):
	# use intermediate folder to avoid file name conflicts
	temp_dir = "JI_test_kmc_database"
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
		subprocess.run(f"mv -n {temp_name} ./{temp_dir}", shell=True )
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
def unwrap_gt_ji_vector(arg):
	return kmc_intersection_count_for_ji(arg[0], arg[1], arg[2])
# get GT JI vector
def get_gt_ji_vector(input_file, other_list, input_k, thread=6):
	pool = Pool(processes=thread)
	Y = np.array(pool.map(unwrap_gt_ji_vector, zip([input_file] * len(other_list), other_list, [input_k] * len(other_list))))
	pool.terminate()
	return Y
# generate GT JI matrix
def gt_ji_matrix(list1, list2, input_k, out_file_name, thread=6):
	lst = []
	row_name = []
	col_name = [os.path.basename(list2[i]) for i in range(len(list2))]
	for i in range(len(list1)):
		name = os.path.basename(list1[i])
		row_name.append(name)
		Y = get_gt_ji_vector(list1[i], list2, input_k, thread=thread)
		lst.append(Y)
	df = pd.DataFrame(lst)
	df.columns = col_name
	df.index = row_name
	df.to_csv(out_file_name, index=True, encoding='utf-8')
	return df



### Create JI_CE object
class JI_CountEstimator(object):
	def __init__(self, n=None, max_prime=9999999999971, ksize=None, input_file_name=None, rev_comp=True,
	             full_kmer=False):
		if n is None:
			raise Exception
		if ksize is None:
			raise Exception
		
		# read parameters
		self.ksize = ksize
		self.maxk = ksize
		self.p = get_prime_lt_x(max_prime)
		self.input_file_name = input_file_name
		self.rev_comp = rev_comp
		self.full_kmer = full_kmer
		
		# initialize sketch of size n with count
		self._mins = [self.p] * n
		self._kmers = [''] * n
		self._all_kmer = dict()
		self._truncated_all_kmer = ['']
		
		# read sequence by "screed" and load them into the sketch (replace parse_file function in L105)
		if self.input_file_name:
			for record in screed.open(self.input_file_name):
				self.add_sequence(record.sequence, self.full_kmer)
	
	def add_sequence(self, seq, update_full=False):
		seq = seq.upper()  # transfer to upper case (just in case)
		seq_split_onlyACTG = re.compile('[^ACTG]').split(seq)  # get rid of non-ATCG letter
		if len(seq_split_onlyACTG) == 1:
			if len(seq) >= self.ksize:
				for kmer in kmers(seq, self.ksize):
					self.add(kmer, update_full)
		else:
			for sub_seq in seq_split_onlyACTG:
				if len(sub_seq) >= self.ksize:  # in case of small chunk
					self.add_sequence(sub_seq, update_full)
	
	def add(self, kmer, update_full=False):
		_mins = self._mins
		_kmers = self._kmers
		# use rev_comp if needed
		if self.rev_comp:
			kmer = min(kmer, khmer.reverse_complement(kmer))
		h = khmer.hash_no_rc_murmur3(kmer)
		# insert into full kmer set
		if update_full:
			_full = self._all_kmer
			if kmer not in _full:
				_full[kmer] = 1
			else:
				_full[kmer] += 1
		# insert into MH sketches reminder of max_prime we use
		h = h % self.p
		# early stop if n sketches are found
		if h >= _mins[-1]:
			return
		# insert kmer into the sketch
		i = bisect.bisect_left(_mins, h)  # find index to insert h
		if _mins[i] == h:  # already in sketch
			return
		else:
			# h not in sketch, insert
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
	
	def ji_vector(self, other_list, thread=6):
		pool = Pool(processes=thread)
		Y = np.array(pool.map(unwrap_jaccard_vector, zip([self] * len(other_list), other_list)))
		pool.terminate()
		return Y
	
	def brute_force_truncation(self, new_ksize):
		if not isinstance(new_ksize, int):
			raise Exception("Input number is not an integer")
		if new_ksize > self.ksize:
			raise Exception("New size must be smaller than %d." % self.ksize)
		elif new_ksize == self.ksize:
			return
		elif new_ksize < self.ksize:
			# data to be updated after the truncation:
			self.ksize = new_ksize
			while self._mins[-1] == self.p:  # rm unused cells, otherwise empty cell (though very rare) has hash value 0
				self._mins.pop()
				self._kmers.pop()
			new_kmers = list(set([x[0:new_ksize] for x in self._kmers]))
			sketch_size = len(new_kmers)
			self._mins = [self.p] * sketch_size
			self._kmers = [''] * sketch_size
			# update
			for i in range(sketch_size):
				self.add(new_kmers[i])  # for MH sketch only
			# clean trailing empty cells in sketches
			while self._mins[-1] == self.p:
				self._mins.pop()
				self._kmers.pop()
			# conditional: truncate the full kmer to current ksize
			if self.full_kmer:
				old_kmers = [x[0:new_ksize] for x in self._all_kmer]
				if self.rev_comp:
					old_kmers = [min(x, khmer.reverse_complement(x)) for x in old_kmers]
				self._truncated_all_kmer = list(set(old_kmers))
			return
	
	def calculate_bias_factor(self, other):
		"""
		Calculate the bias factor from 2 JI_CE object, need to be truncated first
		Will use: 2 truncated full kmer, 2 full kmer, maxk, current ksize
		"""
		if self.ksize != other.ksize:
			raise Exception("different k-mer sizes - cannot compare")
		if self.p != other.p:
			raise Exception("different primes - cannot compare")
		if self.maxk != other.maxk:
			raise Exception("different maxk - cannot compare")
		if not self.full_kmer or not other.full_kmer:
			raise Exception("full kmer not enabled for the CE object")
		
		# use dict to count prefix
		ksmall_intersect = dict()
		for kmer in list(set(self._truncated_all_kmer).intersection(other._truncated_all_kmer)):
			ksmall_intersect[kmer] = 0  # for counting purpose
		ksmall_union = dict()
		for kmer in list(set(self._truncated_all_kmer).union(other._truncated_all_kmer)):
			ksmall_union[kmer] = 0
		
		# count prefix match
		for kmer in list(set(self._all_kmer.keys()).union(other._all_kmer.keys())):
			kmer = kmer[0:self.ksize]  # prefix
			if self.rev_comp:
				kmer = min(kmer, khmer.reverse_complement(kmer))
			if kmer in ksmall_intersect:
				ksmall_intersect[kmer] += 1
				ksmall_union[kmer] += 1
			elif kmer in ksmall_union:
				ksmall_union[kmer] += 1
		
		# bias factor
		numerator = sum(ksmall_intersect.values()) * 1.0 / len(ksmall_intersect)
		denominator = sum(ksmall_union.values()) * 1.0 / len(ksmall_union)
		bias_factor = numerator / denominator
		print(numerator)
		print(denominator)
		return bias_factor
	
	def bias_vector(self, other_list, thread=6):
		pool = Pool(processes=thread)
		Y = np.array(pool.map(unwrap_bias_vector, zip([self] * len(other_list), other_list)))
		pool.terminate()
		return Y


### trunc_JI and est_JI uses same functions to call
# wrapper for trunc_JI from JI_CE object
def unwrap_jaccard_vector(arg):
	return arg[0].est_jaccard(arg[1])
# get JI vector: defined inside the JI_CE object
# get JI matrix
def ji_matrix(list1, list2, out_file_name, thread=6):
	lst = []
	row_name = []
	col_name = [os.path.basename(list2[i].input_file_name) for i in range(len(list2))]
	for i in range(len(list1)):
		name = os.path.basename(list1[i].input_file_name)
		row_name.append(name)
		Y = list1[i].ji_vector(list2, thread=thread)
		lst.append(Y)
	df = pd.DataFrame(lst)
	df.columns = col_name
	df.index = row_name
	df.to_csv(out_file_name, index=True, encoding='utf-8')
	return df



### Bias factor
# wrapper for bias factor from JI_CE object
def unwrap_bias_vector(arg):
	return arg[0].calculate_bias_factor(arg[1])
# bias vector defined inside JI_CE object
# bias matrix
def bias_matrix(list1, list2, out_file_name, thread=6):
	lst = []
	row_name = []
	col_name = [os.path.basename(list2[i].input_file_name) for i in range(len(list2))]
	for i in range(len(list1)):
		name = os.path.basename(list1[i].input_file_name)
		row_name.append(name)
		Y = list1[i].bias_vector(list2, thread=thread)
		lst.append(Y)
	df = pd.DataFrame(lst)
	df.columns = col_name
	df.index = row_name
	df.to_csv(out_file_name, index=True, encoding='utf-8')
	return df


### get JI_CE object lists and array functions
def make_minhash(genome, max_h, prime, ksize, rev_comp, full_kmer):
	MHS = JI_CountEstimator(n=max_h, max_prime=prime, ksize=ksize, input_file_name=genome, rev_comp=rev_comp, full_kmer=full_kmer)
	return MHS

def make_minhash_star(arg):
	return make_minhash(*arg)

def get_ce_lists(input_file, input_k, reverse_comp, full_kmer=False, thread=6, max_h=1000, prime=9999999999971):
	para = Pool(processes=thread)
	genome_sketches = para.map(make_minhash_star, zip(input_file, repeat(max_h), repeat(prime), repeat(input_k), repeat(reverse_comp), repeat(full_kmer)))
	para.close()
	return genome_sketches

def bf_truncation(ce_list, new_ksize):
	for ce in ce_list:
		ce.brute_force_truncation(new_ksize)



### test function
# test ground truth JI
def f1_test_gt_ji():
	subprocess.run(f"echo TTAATTAA > test1.fq", shell=True)
	subprocess.run(f"echo AAATTTAA > test2.fq", shell=True)
	create_kmc_db("test1.fq", 4, thread=6)
	create_kmc_db("test2.fq", 4, thread=6)
	test_ji = kmc_intersection_count_for_ji("test1.fq", "test2.fq", 4)
	os.remove("test1.fq")
	os.remove("test2.fq")
	print("This is for canonical kmers: |A|=3, |B|=4, overlap=2")
	print(test_ji)
	assert test_ji == 2 / 5.0

# test est JI
def f2_test_est_ji():
	E1 = JI_CountEstimator(n=0, ksize=21)
	E2 = JI_CountEstimator(n=0, ksize=21)
	E1._mins = [1, 2, 3, 4, 5]
	E2._mins = [1, 2, 3, 4]
	assert E1.est_jaccard(E2) == 4 / 4.0
	assert E2.est_jaccard(E1) == 4 / 4.0
	E2._mins = [1, 3, 4, 5, 6]
	assert E1.est_jaccard(E2) == 4 / 5.0

# test truncation
def f3_test_truncation():
	E1 = JI_CountEstimator(n=5, ksize=4)
	E2 = JI_CountEstimator(n=5, ksize=4)
	E1.add_sequence("TTAATTAA")
	print(E1._kmers)
	E1.brute_force_truncation(3)
	print(E1._kmers)
	assert E1._kmers == ['AAT', 'TAA']

# test bias factor
def f4_test_bias_calculation():
	# no conflux
	E1 = JI_CountEstimator(n=3, ksize=5, rev_comp=False, full_kmer=True)
	E2 = JI_CountEstimator(n=3, ksize=5, rev_comp=False, full_kmer=True)
	E1.add_sequence("AAAAA", update_full=True)
	E1.add_sequence("TTTTT", update_full=True)
	E1.add_sequence("CCCCC", update_full=True)
	
	E2.add_sequence("AAATT", update_full=True)
	E2.add_sequence("TTTGG", update_full=True)
	E2.add_sequence("CCCGG", update_full=True)
	E1.brute_force_truncation(3)
	E2.brute_force_truncation(3)
	assert E1.calculate_bias_factor(E2) == 1.0  # no dup introduced, get 1
	
	# with dup
	E3 = JI_CountEstimator(n=3, ksize=5, rev_comp=False, full_kmer=True)
	E4 = JI_CountEstimator(n=3, ksize=5, rev_comp=False, full_kmer=True)
	E3.add("AAAAA", update_full=True)
	E3.add("AAATT", update_full=True)
	E3.add("AAAGG", update_full=True)
	
	E4.add("AAACC", update_full=True)
	E4.add("GGGTC", update_full=True)
	E4.add("GGGCC", update_full=True)
	E3.brute_force_truncation(3)
	E4.brute_force_truncation(3)
	assert E3.calculate_bias_factor(E4) == (4 * 1.0) / (6 / 2)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="This script creates training/reference sketches for each FASTA/Q file"
	                                             " listed in the input file.",
	                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-q', '--query', help="Path to the query file.")
	parser.add_argument('-r', '--ref', help="Path to the ref file.")
	parser.add_argument('-f', '--fast', type=str, help="Skipping steps to build GT_TST and CE.h5 file", default="False")
	parser.add_argument('-c', '--rev_comp', type=str, help="Whether to keep the reverse complementary", default="True")
	parser.add_argument('-n', '--num_hashes', type=int, help="Number of hashes to use.", default=1000)
	parser.add_argument('-k', '--k_size', type=int, help="k-mer size", default=60)
	parser.add_argument('-g', '--k_range', type=str, help="k-mer range", default="10-60-5")
	parser.add_argument('-t', '--threads', type=int, help="Number of threads to use",
	                    default=min(64, int(multiprocessing.cpu_count() / 2)))
	parser.add_argument('-p', '--prime', help='Prime (for modding hashes)', default=9999999999971)
	parser.add_argument('-z', '--skip_gt', type=str, help="Skip GT JI calculation", default="False")
	parser.add_argument('-x', '--skip_est', type=str, help="Skip Est JI calculation", default="False")
	parser.add_argument('-y', '--skip_trunc', type=str, help="Skip Trunc JI calculation", default="False")
	parser.add_argument('-o', '--skip_bias', type=str, help="Skip Bias factor calculation", default="False")
	
	args = parser.parse_args()
	num_threads = args.threads
	prime = args.prime  # taking hashes mod this prime
	ksize = args.k_size  # max k used in the analysis
	max_h = args.num_hashes  # number of hashes to use
	input_range = args.k_range  # k range
	k_sizes = parsenumlist(input_range)
	query_file = args.query
	ref_file = args.ref
	rev_comp = args.rev_comp
	fast_mode = args.fast
	skip_gt = args.skip_gt
	skip_est = args.skip_est
	skip_trunc = args.skip_trunc
	skip_bias = args.skip_bias
	
	rev_comp = rev_comp == 'True'
	if rev_comp:
		print("Using true for rev_comp!!!!!!")
	
	fast_mode = fast_mode == 'True'
	if fast_mode:
		print("Skipping all kmc building process and assume they are already built!!!!!!")
	
	skip_gt = skip_gt == 'True'
	if skip_gt:
		print("Skipping all GT JI calculation steps")
	
	skip_est = skip_est == 'True'
	if skip_est:
		print("Skipping all Est JI calculation steps")
	
	skip_trunc = skip_trunc == 'True'
	if skip_trunc:
		print("Skipping all Trunc JI calculation steps")
	
	skip_bias = skip_bias == 'True'
	if skip_bias:
		print("Skipping all Trunc JI calculation steps")
	
	query_file_names = os.path.abspath(query_file)
	if not os.path.exists(query_file_names):
		raise Exception("Input file %s does not exist." % query_file_names)
	ref_file_names = os.path.abspath(ref_file)
	if not os.path.exists(ref_file_names):
		raise Exception("Input file %s does not exist." % query_file_names)
	
	# pipe start
	query_list = check_files(query_file)
	ref_list = check_files(ref_file)
	
	# create kmc database for all files for all k values
	if not fast_mode:
		print("Creating kmc database for JI test")
		if not os.path.exists("JI_test_kmc_database"):
			os.mkdir("JI_test_kmc_database")
		for file in set(query_list).union(set(ref_list)):
			for k in k_sizes:
				print("Running file " + file)
				create_kmc_db(file, k, thread=num_threads)
	
	# run GT_JI by kmc
	if not skip_gt:
		for k in k_sizes:
			out_name = "GroundTruth_JI_k" + str(k) + ".csv"
			gt_ji_matrix(query_list, ref_list, k, out_name, thread=num_threads)
	
	# run trunc_JI by JI_CE object
	if not skip_trunc:
		sketch1 = get_ce_lists(query_list, ksize, rev_comp, thread=num_threads, max_h=max_h, prime=prime)
		sketch2 = get_ce_lists(ref_list, ksize, rev_comp, thread=num_threads, max_h=max_h, prime=prime)
		rev_k_sizes = k_sizes.copy()
		rev_k_sizes.reverse()
		for i in rev_k_sizes:
			bf_truncation(sketch1, i)
			bf_truncation(sketch2, i)
			out_name = "trunc_JI_k" + str(i) + ".csv"
			out1 = ji_matrix(sketch1, sketch2, out_name, thread=num_threads)
	
	# get est_JI
	if not skip_est:
		for i in k_sizes:
			sketch1 = get_ce_lists(query_list, i, rev_comp, thread=num_threads, max_h=max_h, prime=prime)
			sketch2 = get_ce_lists(ref_list, i, rev_comp, thread=num_threads, max_h=max_h, prime=prime)
			out_name = "est_JI_k" + str(i) + ".csv"
			out2 = ji_matrix(sketch1, sketch2, out_name, thread=num_threads)
	
	# get bias factor: same to trunc JI, but separate here for testing purpose
	if not skip_bias:
		sketch1 = get_ce_lists(query_list, ksize, rev_comp, full_kmer=True, thread=num_threads, max_h=max_h, prime=prime)
		sketch2 = get_ce_lists(ref_list, ksize, rev_comp, full_kmer=True, thread=num_threads, max_h=max_h, prime=prime)
		rev_k_sizes = k_sizes.copy()
		rev_k_sizes.reverse()
		if rev_k_sizes[0] == ksize:  # don't need the maxk for bias factor
			rev_k_sizes.remove(ksize)
		for i in rev_k_sizes:
			bf_truncation(sketch1, i)
			bf_truncation(sketch2, i)
			out_name = "bias_factor_k" + str(i) + "_to_k" + str(ksize) + ".csv"
			out1 = bias_matrix(sketch1, sketch2, out_name, thread=num_threads)

################# validate 2 files
#for the 2 files:
# GT_JI ~= Est_JI = 0.676 (no revcompe GT=0.02)
# But Trunc_JI = 0.002, which is the rev_comp=False results of GT, but I do have rev_comp for truncation
#
# os.chdir("./test_JI_bias_factor/validate_trunc_result")
# file1 = 'GCA_001721725.1_ASM172172v1_genomic.fna'
# file2 = 'GCA_001022035.1_ASM102203v1_genomic.fna'
# k_sizes = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
# rev_comp = True
# num_threads=6
# prime=9999999999971
# ksize=60
# max_h=1000
# full_kmer = False
#
# def get_all_kmers(input_file, temp_k, use_rev_comp=True):
# 	temp_dict = dict()
# 	for record in screed.open(input_file):
# 		for kmer in kmers(record.sequence, temp_k):
# 			if use_rev_comp:
# 				kmer = min(kmer, khmer.reverse_complement(kmer))
# 			if kmer in temp_dict:
# 				temp_dict[kmer] += 1
# 			else:
# 				temp_dict[kmer] = 1
# 	return temp_dict
#
# def output_ji(kmer_dict1, kmer_dict2):
# 	intersection = len(set(kmer_dict1.keys()).intersection(set(kmer_dict2.keys())))
# 	union = len(set(kmer_dict1.keys()).union(set(kmer_dict2.keys())))
# 	out_ji = intersection*1.0/union
# 	print(out_ji)
#
# temp_k = 60
# a_kmers = get_all_kmers(file1, temp_k, True)
# b_kmers = get_all_kmers(file2, temp_k, True)
# output_ji(a_kmers, b_kmers)
# # k20 result: 0.6766208404617552
# # k60 result: 0.5005968561188748

# Est_JI at k 20 rev_comp = True
# file_list = [file1, file2]
# temp_k = 60
# sketch1 = get_ce_lists(file_list, temp_k, True, max_h=max_h, prime=prime)
# obj1 = sketch1[0]
# obj2 = sketch1[1]
# obj1.est_jaccard(obj2)
# # k20 result: 0.676 (match)
# # k60 result: 0.492 (match)
#
# #Trunc JI jumps smaller
# temp_k = 60
# sketch1 = get_ce_lists(file_list, temp_k, True, max_h=max_h, prime=prime)
# obj1 = sketch1[0]
# obj2 = sketch1[1]
# obj1.est_jaccard(obj2)
# k60_kmer1 = obj1._kmers
# k60_kmer2 = obj2._kmers
# # est 0.492
# # trun
# bf_truncation(sketch1, 20)
# obj1.est_jaccard(obj2)
# ## trunc k20 result: 0.492 vs 0.676
# k20_kmer1 = obj1._kmers
# k20_kmer2 = obj2._kmers
#
# bf_truncation(sketch1, 10)
# obj1.est_jaccard(obj2)
# k10_kmer1 = obj1._kmers
# k10_kmer2 = obj2._kmers
#
# # check kmers:
# # no shared prefix found: 1 on 1 projection
# len(set(k60_kmer1).intersection(set(k60_kmer2)))
# len(set(k20_kmer1).intersection(set(k20_kmer2)))
# #



