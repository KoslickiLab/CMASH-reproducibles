# usage: python <script> -q <query_file> -r <ref_file>

# import CMash object
from CMash.MinHash import *
# other modules
import multiprocessing
from multiprocessing import Pool
import argparse
import os
import khmer
from itertools import *
import numpy as np
import pandas as pd
import re
from argparse import ArgumentTypeError
import subprocess
import screed
import bisect
import inspect
import copy
import gc


# function defination
#################################################
# workflow related
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


def check_files(input_file):
	"""
	Read all file paths from an input file and check if they exist
	:param input_file:
	:return: a sorted list of file paths in the input file
	"""
	print("Reading paths from %s to a list." %input_file)
	out_list = list()
	with open(input_file, 'r') as temp:
		for line in temp.readlines():
			line = line.strip()
			if not os.path.exists(line):
				raise Exception("Input file %s does not exist." % line)
			out_list.append(os.path.abspath(line))
		out_list = sorted(out_list, key=os.path.basename)
		return (out_list)


#################################################
# CMash Count Estimator (CCE) object definition and adjustment
def estimate_genome_size(genome, ksize):
	"""
	fast kmer cardinality estimate (copied from CMash, but didn't find it from MH)
	:param genome:
	:param ksize:
	:return: estimated kmer cardinality of a genomic file on a given k_size
	"""
	hll = khmer.HLLCounter(0.01, ksize)
	hll.consume_seqfile(genome)
	return hll.estimate_cardinality()


def containment_to_jaccard(C_est, CE1, CE2):
	"""
	transfer CI to JI by their definition
	:param C_est: CI(A, B)
	:param CE1: CE obj1
	:param CE2: CE obj2
	:return: JI(A, B)
	"""
	A = CE1.cardinality
	B = CE2.cardinality
	return C_est * A * 1.0 / (A + B - C_est * A)


# (not a "true" canonical)
def canonical_kmer(kmer):
	"""
	transfer an input kmer into its canomical form
	:param kmer:
	:return:
	"""
	h1 = khmer.hash_no_rc_murmur3(kmer)
	h2 = khmer.hash_no_rc_murmur3(khmer.reverse_complement(kmer))
	if h1 > h2:
		kmer = khmer.reverse_complement(kmer)
	return kmer


def wrapper_est_jaccard_vector(arg):
	"""
	wrapper function for CCE object: est_jaccard function, 1 vs a list
	:param arg:
	:return:
	"""
	return arg[0].est_jaccard(arg[1])


def wrapper_est_jaccard_matrix(list1, list2, out_file_name, threads=8):
	"""
	wrapper function to get a pairwise est_JI matrix for CCE object list1 vs list2
	:param threads:
	:param list1:
	:param list2:
	:param out_file_name:
	:param thread:
	:return: return nothing but an output file with the matrix
	"""
	lst = []
	row_name = []
	col_name = [os.path.basename(list2[i].input_file_name) for i in range(len(list2))]
	for i in range(len(list1)):
		name = os.path.basename(list1[i].input_file_name)
		row_name.append(name)
		Y = list1[i].ji_vector(list2, threads=threads)
		lst.append(Y)
	df = pd.DataFrame(lst)
	df.columns = col_name
	df.index = row_name
	df.to_csv(out_file_name, index=True, encoding='utf-8')
	return


def wrapper_bias_factor_vector(arg):
	"""
	wrapper function for CCE object to calculate bias factor in 1 vs a list manner
	:param arg:
	:return:
	"""
	return arg[0].calculate_bias_factor(arg[1])


def wrapper_bias_factor_matrix(list1, list2, out_file_name, threads=8):
	"""
	wrapper function to get pairwiss bias factor matrix for 2 list of CCE objects
	:param list1:
	:param list2:
	:param out_file_name:
	:param thread:
	:return:
	"""
	lst = []
	row_name = []
	col_name = [os.path.basename(list2[i].input_file_name) for i in range(len(list2))]
	for i in range(len(list1)):
		name = os.path.basename(list1[i].input_file_name)
		row_name.append(name)
		Y = list1[i].bias_vector(list2, threads=threads)
		lst.append(Y)
	df = pd.DataFrame(lst)
	df.columns = col_name
	df.index = row_name
	df.to_csv(out_file_name, index=True, encoding='utf-8')
	return df


def wrapper_streaming_method_ji_from_ci_vector(arg):
	"""
	wrapper function for streaming method CI -> JI in CCE obj to handle 1 vs a list of objects
	:param arg:
	:return:
	"""
	return arg[0].stream_cmash_for_ji(arg[1])


def wrapper_streaming_method_ji_from_ci_matrix(list1, list2, out_file_name, threads=8):
	"""
	wrapper functino to generate a JI matrix for 2 lists of CCE objects
	:param list1:
	:param list2:
	:param out_file_name:
	:param threads:
	:return:
	"""
	lst = []
	row_name = []
	col_name = [os.path.basename(list2[i].input_file_name) for i in range(len(list2))]
	for i in range(len(list1)):
		name = os.path.basename(list1[i].input_file_name)
		row_name.append(name)
		Y = list1[i].cmash_stream_ji_vector(list2, threads=threads)
		lst.append(Y)
	df = pd.DataFrame(lst)
	df.columns = col_name
	df.index = row_name
	df.to_csv(out_file_name, index=True, encoding='utf-8')
	return df


def wrapper_ground_truth_ji_vector(arg):
	"""
	wrapper function for ground truth JI calculation for CCE objects 1 vs a list
	:param arg:
	:return:
	"""
	return arg[0].ground_truth_jaccard(arg[1])


def wrapper_ground_truth_ji_matrix(list1, list2, out_file_name, threads=8):
	"""
	wrapper functino to generate a GT_JI matrix for 2 lists of CCE objects
	:param list1:
	:param list2:
	:param out_file_name:
	:param threads:
	:return:
	"""
	lst = []
	row_name = []
	col_name = [os.path.basename(list2[i].input_file_name) for i in range(len(list2))]
	for i in range(len(list1)):
		name = os.path.basename(list1[i].input_file_name)
		row_name.append(name)
		Y = list1[i].ground_truth_ji_vector(list2, threads=threads)
		lst.append(Y)
	df = pd.DataFrame(lst)
	df.columns = col_name
	df.index = row_name
	df.to_csv(out_file_name, index=True, encoding='utf-8')
	return

# the main CCE object (inherit from CE)
class CMash_CountEstimator(CountEstimator):
	"""
	CE object for CMash manuscrip
	Inheriate CountEstimator object from CMash.MinHash
	"""
	
	
	# initialization
	def __init__(self, n=None, max_prime=9999999999971., ksize=None, input_file_name=None, rev_comp=True, full_kmer=False, threads=8):
		CountEstimator.__init__(self, n=n, max_prime=max_prime, ksize=ksize, input_file_name=input_file_name, rev_comp=True, save_kmers='n', hash_list=None)
		# parameters defined in CE: self.p, self.ksize, self.input_file_name,
		
		if n is None:
			raise Exception
		if ksize is None:
			raise Exception
		
		# read parameters
		self.maxk = ksize   # max k value for the truncation
		self.rev_comp = rev_comp # binary indicator: use canonical k-mer, default Yes
		self.full_kmer = full_kmer # binary indicator: store full k-mer sets besides the n-sketch for CMash validation
		self.threads=threads
		
		# initialize kmer/hash vector (set save_kmer != y in __init__, so do it manually here)
		self._mins = [self.p] * n   # hash sketch
		self._kmers = [''] * n  # kmer sketch
		self._all_kmer = dict() # for bias factor calculation, store full k-mer set of a given file
		self._truncated_all_kmer = [''] # for bias factor calculation, store truncated full k-mer set of a given file
		
		# build bottom sketch
		if self.input_file_name:
			self.cardinality = estimate_genome_size(input_file_name, self.ksize)
			for record in screed.open(self.input_file_name):
				self.add_sequence(record.sequence, rev_comp=self.rev_comp)
	
	
	# function: add_full_kmer_set
	### single add operation
	def add_to_full(self, kmer, rev_comp=False):
		# similar to "add" function, add a kmer into full kmer set
		_full = self._all_kmer
		# reverse comp
		if rev_comp:
			kmer=canonical_kmer(kmer)
		# insert into dict
		if kmer not in _full:
			_full[kmer] = 1
		else:
			_full[kmer] += 1
	### main function
	### this step can be merged into "build bottom sketch" section, but might be more convenient in separate cases when we don't need full sets
	def add_full_kmer_set(self):
		# change the indicator
		self.full_kmer=True
		
		for record in screed.open(self.input_file_name):
			# shadow add_sequence function, can be merged together though
			seq = record.sequence
			seq = seq.upper()
			seq_split_onlyACTG = notACTG.split(seq)
			for sub_seq in seq_split_onlyACTG:
				if sub_seq:
					for kmer in kmers(sub_seq, ksize=self.ksize):
						self.add_to_full(kmer, rev_comp=self.rev_comp)
		
		# it's not necessary to keep the dict "self._all_kmers" by now, but I still keep it in case we need a weighted count in the future
		self._truncated_all_kmer = list(self._all_kmer.keys())
	
	
	# estimate JI between 2 CE objects, slightly dif from CMash (same idea)
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
			processed += 1
			if mins1[i] < mins2[j]:
				i += 1
			elif mins1[i] > mins2[j]:
				j += 1
			elif mins1[i] == mins2[j]:
				# a match from A overlap B
				overlap += 1
				i += 1
				j += 1
		print(overlap)
		print(processed)
		est_ji = overlap * 1.0 / processed
		return est_ji
	
	
	# generate a JI vector for self vs a list of CCE objects, "unwrap_jaccard_vector" from MH
	def ji_vector(self, other_list, threads=8):
		pool = Pool(processes=threads)
		Y = np.array(pool.map(wrapper_est_jaccard_vector, zip([self] * len(other_list), other_list)))
		pool.terminate()
		return Y
	
	
	# brute-forcely truncate k-mer sketches of self, this step is NOT reversible!!! (can be, but not necessary)
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
			if self.input_file_name:
				self.cardinality = estimate_genome_size(self.input_file_name, new_ksize)
			while self._mins[-1] == self.p:  # rm unused cells, otherwise empty cell (though very rare) has hash value 0
				self._mins.pop()
				self._kmers.pop()
			new_kmers = list(set([x[0:new_ksize] for x in self._kmers]))
			sketch_size = len(new_kmers)
			self._mins = [self.p] * sketch_size
			self._kmers = [''] * sketch_size
			# update
			for i in range(sketch_size):
				self.add(new_kmers[i], rev_comp=self.rev_comp)  # for MH sketch only
			# clean trailing empty cells in sketches
			while self._mins[-1] == self.p:
				self._mins.pop()
				self._kmers.pop()
			# conditional: truncate the full kmer to current ksize
			if self.full_kmer:
				trunc_kmers = [x[0:new_ksize] for x in self._all_kmer]
				if self.rev_comp:
					trunc_kmers = [canonical_kmer(x) for x in trunc_kmers]
				self._truncated_all_kmer = list(set(trunc_kmers))
			return
	
	
	# calculate bias factor
	def calculate_bias_factor(self, other):
		"""
		Calculate the bias factor for 2 CCE object, need to be truncated first
		Input: 2 truncated full kmer, 2 full kmer, maxk, current ksize
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
				kmer = canonical_kmer(kmer)
			if kmer in ksmall_intersect:
				ksmall_intersect[kmer] += 1
				ksmall_union[kmer] += 1
			elif kmer in ksmall_union:
				ksmall_union[kmer] += 1
		
		# bias factor
		if len(ksmall_intersect) == 0:
			numerator = 0
		else:
			numerator = sum(ksmall_intersect.values()) * 1.0 / len(ksmall_intersect)
		denominator = sum(ksmall_union.values()) * 1.0 / len(ksmall_union)
		bias_factor = numerator / denominator
		print(numerator)
		print(denominator)
		return bias_factor
		

	# generate a bias factor vector for self vs a list of CCE objects
	def bias_vector(self, other_list, threads=8):
		pool = Pool(processes=threads)
		Y = np.array(pool.map(wrapper_bias_factor_vector, zip([self] * len(other_list), other_list)))
		pool.terminate()
		return Y


	# streaming method to est JI from CI
	def stream_cmash_for_ji(self, other):
		"""
		Est JI from CI by streaming method
		We can only get truncated_JI in this way. It's not correct to directly truncate the hash lists from 2 CCE objects (wrong sampling method)
		:param other:
		:return:
		"""
		if self.ksize != other.ksize:
			raise Exception("different k-mer sizes - cannot compare")
		if self.p != other.p:
			raise Exception("different primes - cannot compare")
		
		A_kmers = self._kmers
		A_matches = dict()
		for kmer in A_kmers:
			if self.rev_comp:
				kmer = canonical_kmer(kmer)
			A_matches[kmer] = 0  # count purpose
		
		# streaming all other kmers for CI (can't do JI directly)
		for record in screed.open(other.input_file_name):
			seq = record.sequence
			seq = seq.upper()
			seq_split_onlyACTG = re.compile('[^ACTG]').split(seq)
			for sub_seq in seq_split_onlyACTG:
				for i in range(len(sub_seq) - self.ksize + 1):
					# enumerate all kmers
					kmer = sub_seq[i:i + self.ksize]
					if self.rev_comp:
						kmer = canonical_kmer(kmer)
					if kmer in A_matches:
						A_matches[kmer] = 1
		
		# return results
		C_est = np.sum(list(A_matches.values())) / len(A_kmers)
		J_est = containment_to_jaccard(C_est, self, other)
		print(C_est)
		print(J_est)
		return J_est
	
	
	# generate a trunc_JI vector (from CI by streaming method) for self vs a list of CCE objects
	def cmash_stream_ji_vector(self, other_list, threads=8):
		pool = Pool(processes=threads)
		Y = np.array(pool.map(wrapper_streaming_method_ji_from_ci_vector, zip([self] * len(other_list), other_list)))
		pool.terminate()
		return Y


	# ground truth JI: discard kmc which takes toooooo many storage spaces
	def ground_truth_jaccard(self, other):
		if self.ksize != other.ksize:
			raise Exception("different k-mer sizes - cannot compare")
		if self.p != other.p:
			raise Exception("different primes - cannot compare")
		if not self.full_kmer or not other.full_kmer:
			raise Exception("full kmer not enabled for the CE object")
		
		_inter = set(self._truncated_all_kmer).intersection(other._truncated_all_kmer)
		_union = set(self._truncated_all_kmer).union(other._truncated_all_kmer)
		print(len(_inter))
		print(len(_union))
		return len(_inter)*1.0/len(_union)
	
	
	# get a GT_JI vector
	def ground_truth_ji_vector(self, other_list, threads=8):
		pool = Pool(processes=threads)
		Y = np.array(pool.map(wrapper_ground_truth_ji_vector, zip([self] * len(other_list), other_list)))
		pool.terminate()
		return Y


#################################################
# CMash analysis workflow
def make_cce_object(genome, max_h, prime, ksize, rev_comp, threads):
	"""
	a wrapper function to create CCE object
	:param genome:
	:param max_h:
	:param prime:
	:param ksize:
	:param rev_comp:
	:param threads:
	:return:
	"""
	CCE = CMash_CountEstimator(n=max_h, max_prime=prime, ksize=ksize, input_file_name=genome, rev_comp=rev_comp, threads=threads)
	return CCE


def wrapper_make_cce_star(arg):
	"""
	a wrapper function to call "make_cce_object" by batch
	:param arg:
	:return:
	"""
	return make_cce_object(*arg)


def get_cce_lists(input_file, input_k, rev_comp, threads=8, max_h=1000, prime=9999999999971):
	"""
	build a list of CCE objects from an input file
	:param input_file:
	:param input_k:
	:param rev_comp:
	:param threads:
	:param max_h:
	:param prime:
	:return:
	"""
	para = Pool(processes=threads)
	genome_sketches = para.map(wrapper_make_cce_star, zip(input_file, repeat(max_h), repeat(prime), repeat(input_k), repeat(rev_comp), repeat(threads)))
	para.close()
	return genome_sketches


def cce_bf_trunc(cce_obj, new_ksize):
	"""
	wrapper function to do bf_trunc of a given CCE object
	:param new_ksize:
	:param cce_obj:
	:return:
	"""
	cce_obj.brute_force_truncation(new_ksize)
	return cce_obj


def cce_bf_trunc_star(arg):
	"""
	just a wrapper to pass 2 para to cce_bf_trunc function
	:param arg:
	:return:
	"""
	return cce_bf_trunc(*arg)


def cce_load_full_kmer(cce_obj):
	"""
	wrapper function for a CCE to call load full kmer function
	This will create a copy of CCE obj
	:param cce_obj:
	:return:
	"""
	cce_obj.add_full_kmer_set()
	return cce_obj

	
#################################################
# test fucntion: validate pipeline performance on known data

# test GT_JI
def f1_test_gt_ji():
	subprocess.run(f"echo '>test_seq\nTTAATTAA' > test1.fq", shell=True)
	subprocess.run(f"echo '>test_seq\nAAATTTAA' > test2.fq", shell=True)
	obj1 = make_cce_object(genome="test1.fq", max_h=10, prime=9999999999971, ksize=4, rev_comp=True, threads=1)
	obj1.add_full_kmer_set()
	obj2 = make_cce_object(genome="test2.fq", max_h=10, prime=9999999999971, ksize=4, rev_comp=True, threads=1)
	obj2.add_full_kmer_set()
	os.remove("test1.fq")
	os.remove("test2.fq")
	#print("This is for canonical kmers: |A|=3, |B|=4, overlap=2")
	temp_gt_ji = obj1.ground_truth_jaccard(obj2)
	assert temp_gt_ji == 2 / 5.0


# test est_JI
def f2_test_est_ji():
	E1 = CMash_CountEstimator(n=0, ksize=5)
	E2 = CMash_CountEstimator(n=0, ksize=5)
	E1._mins = [1, 2, 3, 4, 5]
	E2._mins = [1, 2, 3, 4]
	assert E1.est_jaccard(E2) == 4 / 4.0    #sample until min(|A|, |B|), which is 4 here
	assert E2.est_jaccard(E1) == 4 / 4.0
	E2._mins = [1, 3, 4, 5, 6]
	assert E1.est_jaccard(E2) == 4 / 5.0


# test truncation
def f3_test_truncation():
	E1 = CMash_CountEstimator(n=4, ksize=4)
	E1.add_sequence("TTAATTAA")
	print(E1._kmers)
	E1.brute_force_truncation(3)
	print(E1._kmers)
	assert E1._kmers == ['TTA', 'ATT']  #note: current canonical is based on min(h1, h2), otherwise is ['ATT', 'TAA']


# test bias factor
def f4_test_bias_calculation():
	# no conflux
	E1 = CMash_CountEstimator(n=3, ksize=5, rev_comp=False, full_kmer=True)
	E2 = CMash_CountEstimator(n=3, ksize=5, rev_comp=False, full_kmer=True)
	# _mins and _kmers won't be used, just to make object function viable
	E1._mins = [1,2,3]
	E1._kmers = ["AAAAA", "TTTTT", "CCCCC"]
	E2._mins = [1, 2, 3]
	E2._kmers = ["AAATT", "TTTGG", "CCCGG"]
	#####################
	E1._all_kmer = {k:v for k, v in zip(["AAAAA", "TTTTT", "CCCCC"], [1,1,1])}
	E2._all_kmer = {k:v for k, v in zip(["AAATT", "TTTGG", "CCCGG"], [1,1,1])}
	E1.brute_force_truncation(3)
	E2.brute_force_truncation(3)
	assert E1.calculate_bias_factor(E2) == 1.0  # no dup introduced, get 1
	
	# with dup
	E1 = CMash_CountEstimator(n=3, ksize=5, rev_comp=False, full_kmer=True)
	E2 = CMash_CountEstimator(n=3, ksize=5, rev_comp=False, full_kmer=True)
	# _mins and _kmers won't be used, just to make object function viable
	E1._mins = [1, 2, 3]
	E1._kmers = ["AAAAA", "TTTTT", "CCCCC"]
	E2._mins = [1, 2, 3]
	E2._kmers = ["AAATT", "TTTGG", "CCCGG"]
	#####################
	E1._all_kmer = {k: v for k, v in zip(["AAAAA", "AAATT", "AAAGG"], [1, 1, 1])}
	E2._all_kmer = {k: v for k, v in zip(["AAACC", "GGGTC", "GGGCC"], [1, 1, 1])}
	E1.brute_force_truncation(3)
	E2.brute_force_truncation(3)
	assert E1.calculate_bias_factor(E2) == (4 * 1.0) / (6 / 2)



#################################################
# CMash analysis

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="CMash reproducible analysis",
	                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-q', '--query', help="Path to the query file.")
	parser.add_argument('-r', '--ref', help="Path to the ref file.")
	parser.add_argument('-c', '--rev_comp', type=str, help="Whether to keep the reverse complementary", default="True")
	parser.add_argument('-n', '--num_hashes', type=int, help="Number of hashes to use.", default=2000)
	parser.add_argument('-k', '--k_size', type=int, help="k-mer size", default=60)
	parser.add_argument('-g', '--k_range', type=str, help="k-mer range", default="15-60-5")
	parser.add_argument('-t', '--threads', type=int, help="Number of threads to use",
	                    default=min(24, int(multiprocessing.cpu_count() / 2)))
	parser.add_argument('-p', '--prime', help='Prime (for modding hashes)', default=9999999999971)
	# for test purpose
	parser.add_argument('-z', '--skip_gt', type=str, help="Skip GT JI calculation", default="False")
	parser.add_argument('-x', '--skip_est', type=str, help="Skip Est JI calculation", default="False")
	parser.add_argument('-y', '--skip_trunc', type=str, help="Skip Trunc JI calculation", default="False")
	parser.add_argument('-o', '--skip_bias', type=str, help="Skip Bias factor calculation", default="False")
	
	
	
	# read parameters
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
	skip_gt = args.skip_gt
	skip_est = args.skip_est
	skip_trunc = args.skip_trunc
	skip_bias = args.skip_bias
	
	rev_comp = rev_comp == 'True'
	if rev_comp:
		print("Using true for rev_comp!!!!!!")
	
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
		print("Skipping all bias factor calculation steps")
	
	query_file_names = os.path.abspath(query_file)
	if not os.path.exists(query_file_names):
		raise Exception("Input file %s does not exist." % query_file_names)
	ref_file_names = os.path.abspath(ref_file)
	if not os.path.exists(ref_file_names):
		raise Exception("Input file %s does not exist." % query_file_names)
	
	
	
	# pipe start
	### list of files to build CCE objects
	query_list = check_files(query_file)
	ref_list = check_files(ref_file)
	
	
	### Est_JI and GT_JI can be done together:
	if not skip_est or not skip_gt:
		for i in k_sizes:
			# build CCE obj lists
			sketch1 = get_cce_lists(input_file=query_list, input_k=i, rev_comp=rev_comp, threads=num_threads,
			                        max_h=max_h, prime=prime)
			sketch2 = get_cce_lists(input_file=ref_list, input_k=i, rev_comp=rev_comp, threads=num_threads,
			                        max_h=max_h, prime=prime)
			
			# est_JI
			if not skip_est:
				print("Calculating est_JI for k "+str(i))
				out_name = "est_JI_k" + str(i) + ".csv"
				wrapper_est_jaccard_matrix(sketch1, sketch2, out_name, threads=num_threads)
				
			# GT_JI
			if not skip_gt:
				print("Calculating GT_JI for k "+str(i))
				out_name = "GroundTruth_JI_k" + str(i) + ".csv"
				# read full kmers into the CCE objects
				with Pool(num_threads) as P:
					full_sketch1 = P.map(cce_load_full_kmer, sketch1)
					full_sketch2 = P.map(cce_load_full_kmer, sketch2)
				# get GT
				wrapper_ground_truth_ji_matrix(full_sketch1, full_sketch2, out_name, threads=num_threads)
				del full_sketch1, full_sketch2
				
			# delete obj to release MEM
			del sketch1, sketch2
			gc.collect()
				
				
	### Trunc_JI and bias_factor can be done together:
	if not skip_bias or not skip_trunc:
		# build sketch for maxk
		sketch1 = get_cce_lists(input_file=query_list, input_k=ksize, rev_comp=rev_comp, threads=num_threads,
		                        max_h=max_h, prime=prime)
		sketch2 = get_cce_lists(input_file=ref_list, input_k=ksize, rev_comp=rev_comp, threads=num_threads,
		                        max_h=max_h, prime=prime)
		# load full kmer set
		with Pool(num_threads) as P:
			full_sketch1 = P.map(cce_load_full_kmer, sketch1)
			full_sketch2 = P.map(cce_load_full_kmer, sketch2)
		# release MEM
		del sketch1, sketch2
		gc.collect()
		# reverse ksize: trunc up -> down
		rev_k_sizes = k_sizes.copy()
		rev_k_sizes.reverse()
		if rev_k_sizes[0] == ksize:
			# no truncation
			rev_k_sizes.remove(ksize)
		for i in rev_k_sizes:
			with Pool(num_threads) as P:
				trunc_full_sketch1 = P.map(cce_bf_trunc_star, zip(full_sketch1, repeat(i)))
				trunc_full_sketch2 = P.map(cce_bf_trunc_star, zip(full_sketch2, repeat(i)))
			if not skip_trunc:
				print("Calculating trunc_JI for k "+str(i))
				out_name = "trunc_JI_k" + str(i) + ".csv"
				wrapper_streaming_method_ji_from_ci_matrix(trunc_full_sketch1, trunc_full_sketch2, out_name, threads=num_threads)
			if not skip_bias:
				print("Calculatin bias factor for k "+str(i))
				out_name = "bias_factor_k" + str(i) + "_to_k" + str(ksize) + ".csv"
				wrapper_bias_factor_matrix(trunc_full_sketch1, trunc_full_sketch2, out_name, threads=num_threads)
			# release MEM
			del trunc_full_sketch1, trunc_full_sketch2
			gc.collect()

			


