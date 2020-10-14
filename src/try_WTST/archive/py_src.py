import time
import os
import screed
import re
import bisect
import khmer
import h5py
import numpy as np
import copy
import multiprocessing
from multiprocessing import Pool
from argparse import ArgumentTypeError
import subprocess

# test import
def test_run():
	print("Happy coding : )")


# time a function
def run_time(func):
	start_time = time.time()
	eval(func)
	end_time = time.time()
	dif_time = end_time - start_time
	print("Running time is %d" %dif_time)


# check paths from input files and store them in a list
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


# prime checker: copied from: https://github.com/dkoslicki/CMash/blob/1485a0cde4685c52ae8cfb8bdaa5d1bad77eaac3/CMash/MinHash.py#L821
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
# find a prime lt input value
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


# generate kmers for a sequence
def kmers(seq, ksize):
	for i in range(len(seq) - ksize + 1):
		yield seq[i:i+ksize]


# calculate weighted JI from 2 hash lists and 2 count (weight) lists
def weighted_jaccard(mins1, mins2, counts1, counts2):
	"""
	Returns weighted jaccard index:
	J_{W}=\frac{\sum_{d\subseteq D}min(W_A(d), W_B(d))}{\sum_{d\subseteq D}max(W_A(d), W_B(d))}
	:param mins1: minhash value for CE object1 in a list
	:param mins2: minhash value for CE object2 in a list
	:param counts1: weight for min1
	:param counts2: weight for min2
	:return:
	"""
	sum_min = 0 #for numerator
	sum_max = 0 #for denominator
	i = 0
	j = 0
	processed = 0
	# if a MH value appears only in one list, then min(W(A), W(B)) = 0
	# so we only count the other one if there is no match
	while processed < min(len(mins1), len(mins2)):  # loop stop when min(|A|, |B|) sample size reached
		if mins1[i] < mins2[j]:
			sum_max += counts1[i] #W(B)=0, so add W(A)
			i += 1
			processed += 1
		elif mins1[i] > mins2[j]:
			sum_max += counts2[j] #W(A)=0, so add W(B)
			j += 1
			processed += 1
		elif mins1[i] == mins2[j]:
			#skip the value = max_prime check becaue their weights are all 0, so min/max are 0s
			sum_min += min(counts1[i], counts2[j])
			sum_max += max(counts1[i], counts2[j])
			i += 1
			j += 1
			processed += 1
	WJI = sum_min / float(sum_max)
	print(sum_min) #just for checking purpose
	print(sum_max)
	return WJI
# parallel wrapper for WJI calculation
def unwrap_jaccard_vector(arg):
	return arg[0].est_weighted_jaccard(arg[1])


# parse k range like 10-60-5
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


# bottom n-sketch MinHash implementation
class CountEstimator(object):
	"""
	n: number of sketches to keep
	max_prime: pick a prime <= this nummber
	ksize: kmer length
	input_file: input fasta file / fastq
	hash_list: if a pre-defined hash list is required
	rev_comp: use the min(seq, ser_reverse)
	"""
	def __init__(self, n=None, max_prime=9999999999971, ksize=None, input_file_name=None, hash_list=None, rev_comp=True):
		if n is None:
			raise Exception
		if ksize is None:
			raise Exception

		# read parameters
		self.ksize = ksize
		self.hash_list = hash_list
		self.p = get_prime_lt_x(max_prime)
		self.input_file_name = input_file_name
		self.rev_comp = rev_comp

		# initialize sketch of size n with count
		self._mins = [self.p] * n
		self._counts = [0] * n
		self._kmers = [''] * n
		self._true_num_kmers = 0

		# read sequence by "screed" and load them into the sketch (replace parse_file function in L105)
		if self.input_file_name:
			for record in screed.open(self.input_file_name):
				self.add_sequence(record.sequence, self.rev_comp)

		# at here, this object will read all sequences from the file and store them in the minhash sketch with counts

	# add sequence into sketch
	def add_sequence(self, seq, rev_comp):
		seq = seq.upper() #transfer to upper case (just in case)
		seq_split_onlyACTG = re.compile('[^ACTG]').split(seq) # get rid of non-ATCG letter
		if len(seq_split_onlyACTG) == 1:
			if len(seq) >= self.ksize:
				for kmer in kmers(seq, self.ksize):
					self.add(kmer, 1, rev_comp)
		else:
			for sub_seq in seq_split_onlyACTG:
				if len(sub_seq)>=self.ksize:        #in case of small chunk
					self.add_sequence(sub_seq, rev_comp)

	# add a kmer into the sketch
	def add(self, kmer, weight, rev_comp):
		_mins = self._mins
		_counts = self._counts
		_kmers = self._kmers
		# use rev_comp if needed
		if rev_comp:
			h1 = khmer.hash_no_rc_murmur3(kmer)
			h2 = khmer.hash_no_rc_murmur3(khmer.reverse_complement(kmer))
			h = min(h1, h2)
			if h == h2:
				kmer = khmer.reverse_complement(kmer)
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
			_counts[i] += weight
		else:
			#h not in sketch, insert
			_mins.insert(i, h)
			_counts.insert(i, weight)
			_kmers.insert(i, kmer)
			_mins.pop()
			_counts.pop()
			_kmers.pop()
			return

	# calculate est weighted Jaccard index
	def est_weighted_jaccard(self, other):
		if self.ksize != other.ksize:
			raise Exception("different k-mer sizes - cannot compare")
		if self.p != other.p:
			raise Exception("different primes - cannot compare")
		est_wji = weighted_jaccard(self._mins, other._mins, self._counts, other._counts)
		return est_wji

	# export object using hdf5
	def export(self, export_file_name):
		fid = h5py.File(export_file_name, 'w')
		grp = fid.create_group("CountEstimator")
		mins_data = grp.create_dataset("mins", data=self._mins)
		counts_data = grp.create_dataset("counts", data=self._counts)
		kmer_data = grp.create_dataset("kmers", data=[np.string_(kmer) for kmer in self._kmers])
		grp.attrs['class'] = np.string_("CountEstimator")
		grp.attrs['filename'] = np.string_(self.input_file_name)
		grp.attrs['ksize'] = self.ksize
		grp.attrs['prime'] = self.p
		grp.attrs['true_num_kmers'] = self._true_num_kmers
		fid.close()

	# jaccard vector for many comparisons
	def wji_vector(self, other_list, threads=4):
		pool = Pool(processes=threads)
		Y = np.array(pool.map(unwrap_jaccard_vector, zip([self] * len(other_list), other_list)))
		pool.terminate()
		return Y

	# Non-WTST truncation: directly truncate the kmer list in the object
	def brute_force_truncation(self, new_ksize):
		"""
		Directly truncate the kmer list from a CE object, note this step is NOT reversible!!!
		No WTST involved, this step is only to test the correctness of WTST truncation.
		:param new_ksize: a smaller k than self.ksize
		:return: updated CE object
		"""
		if not isinstance(new_ksize, int):
			raise Exception("Input number is not an integer")
		if new_ksize > self.ksize:
			raise Exception("New size must be smaller than %d." %self.ksize)
		elif new_ksize == self.ksize:
			return
		elif new_ksize < self.ksize:
			# data to be updated after the truncation:
			self.ksize = new_ksize
			new_kmers = [x[0:new_ksize] for x in self._kmers]
			old_counts = self._counts.copy()
			sketch_size = len(new_kmers)
			self._mins = [self.p] * sketch_size
			self._counts = [0] * sketch_size
			self._kmers = [''] * sketch_size
			# update
			for i in range(sketch_size):
				self.add(new_kmers[i], old_counts[i], self.rev_comp)
			while self._mins[-1] == self.p: # rm unused cells
				self._mins.pop()
				self._counts.pop()
				self._kmers.pop()
			return


# export multiple CEs to single hdf5
def export_multiple_to_single_hdf5(CEs, export_file_name):
	"""
	This will take a list of count estimators and export them to a single, large HDF5 file
	:param CEs: list of Count Estimators
	:param file_name: output file name
	:return: None
	"""
	fid = h5py.File(export_file_name, 'w')
	grp = fid.create_group("CountEstimators")
	for CE in CEs:
		try:
			subgrp = grp.create_group(os.path.basename(CE.input_file_name))  # the key of a subgroup is the basename (not the whole file)
			mins_data = subgrp.create_dataset("mins", data=CE._mins)
			counts_data = subgrp.create_dataset("counts", data=CE._counts)
			if CE._kmers is not None:
				kmer_data = subgrp.create_dataset("kmers", data=[np.string_(kmer) for kmer in CE._kmers])
			subgrp.attrs['class'] = np.string_("CountEstimator") # use "_" to apply to various data types
			subgrp.attrs['filename'] = np.string_(CE.input_file_name)  # But keep the full file name on hand
			subgrp.attrs['ksize'] = CE.ksize
			subgrp.attrs['prime'] = CE.p
			subgrp.attrs['true_num_kmers'] = CE._true_num_kmers
		except ValueError:
			fid.close()
			raise Exception("It appears that the training file name %s exists twice in the input data. Please make sure all names are unique (i.e. remove duplicates) and try again." % CE.input_file_name)
	fid.close()
# import multiple CEs from a single hdf5
def import_multiple_from_single_hdf5(file_name, import_list=None):
	"""
	This function will import multiple count estimators stored in a single HDF5 file.
	:param file_name: file name for the single HDF5 file
	:param import_list: List of names of files to import
	:return: a list of Count Estimators
	"""
	CEs = list()
	fid = h5py.File(file_name, 'r')
	if "CountEstimators" not in fid:
		fid.close()
		raise Exception("This function imports a single HDF5 file containing multiple sketches."
						" It appears you've used it on a file containing a single sketch."
						"Try using import_single_hdf5 instead")

	grp = fid["CountEstimators"]
	if import_list:
		iterator = [os.path.basename(item) for item in import_list]
	else:
		iterator = grp.keys()

	iterator = sorted(iterator, key=os.path.basename)  # sort so that we know the order of the input

	for key in iterator:
		if key not in grp:
			fid.close()
			raise Exception("The key " + key + " is not in " + file_name)

		subgrp = grp[key]
		file_name = subgrp.attrs['filename']
		ksize = subgrp.attrs['ksize']
		prime = subgrp.attrs['prime']
		mins = subgrp["mins"][...]
		counts = subgrp["counts"][...]
		true_num_kmers = subgrp.attrs["true_num_kmers"]
		CE = CountEstimator(n=len(mins), max_prime=3, ksize=ksize) #create an empty one and then attach data in
		CE.p = prime
		CE._mins = mins
		CE._counts = counts
		CE._true_num_kmers = true_num_kmers
		CE.input_file_name = file_name
		if "kmers" in subgrp:
			temp_kmers = subgrp["kmers"][...]
			CE._kmers = [kmer.decode('utf-8') for kmer in temp_kmers]
		else:
			CE._kmers = None

		CEs.append(CE)

	fid.close()
	return(CEs)


# helper function for batch construction
class CE_map(object):
	def __init__(self, n, max_prime, ksize):
		self.n = n
		self.max_prime = max_prime
		self.ksize = ksize
	def __call__(self, input_file):
		return CountEstimator(n=self.n, max_prime=self.max_prime, ksize=self.ksize, input_file_name=input_file)
# batch construction for CE objects
def compute_multiple(n=None, max_prime=9999999999971., ksize=None, input_files_list=None, num_threads=1):
	"""
	Parallelly create CEs for a list of input files
	:param n: hash numbers
	:param max_prime:
	:param ksize:
	:param input_files_list:
	:param num_threads:
	:return:
	"""
	if n is None:
		raise Exception
	if ksize is None:
		raise Exception
	if input_files_list is None:
		raise Exception
	pool = Pool(processes=num_threads)
	CEs = pool.map(CE_map(n, max_prime, ksize), input_files_list)
	pool.close()
	return CEs


### calculate ground truth WJI for
# kmc enumerate all kmers
def kmc_count(input_file_name, output_file_name, kmer_size: int, threads=1, rev_comp=False):
	"""
	Call kmc to get all kmers in an input file and store them for future usage
	:param input_file: a fasta file or bam file
	:param out_file_name:
	:return: write down files containing those kmers
	"""
	input_types = ['-fm', '-fq', '-fa', '-fbam']
	success = False
	use_canonical="-b"
	if rev_comp:
		use_canonical=" "
	for input_type in input_types:
		res = subprocess.run(
			f"kmc -k{kmer_size} {input_type} {use_canonical} -r -t{threads} -ci0 -j{output_file_name}.log {input_file_name} {output_file_name} ."
			, shell=True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
		if res.returncode == 0:
			success = True
			break
	if not success:
		raise Exception(
			f"Unknown sequence format: must be one of multifasta, fastq, fasta, or BAM (gzipped or uncompressed). Culprit file is {input_file_name}. Command was {res.args}")
# calculate ground truth WJI
def kmc_intersection_count_for_wji(kmc_input_file1: str, kmc_input_file2: str, keep_intermediate=False):
	"""
	Count sum of min for overlaps and max for all kmers for WJI
	"""

	# numerator: intersect -ocmin
	res = subprocess.run(f"kmc_tools simple {kmc_input_file1} -ci1 {kmc_input_file2} -ci1 intersect temp_intersect_min -ocmin",
						shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	if res.returncode != 0:
		raise Exception(f"The command {res.args} failed to run and returned the returncode={res.returncode}")
	res = subprocess.run(f"kmc_dump temp_intersect_min temp_dump_min; cat temp_dump_min | cut -f 2 | paste -sd+ - | bc", shell=True,
						 capture_output=True)
	if res.returncode != 0:
		raise Exception(f"The command {res.args} failed to run and returned the returncode={res.returncode}")
	wji_numerator = int(res.stdout)
	print(wji_numerator)

	# denominator: union -ocmax
	res = subprocess.run(
		f"kmc_tools simple {kmc_input_file1} -ci1 {kmc_input_file2} -ci1 union temp_union_max -ocmax",
		shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	if res.returncode != 0:
		raise Exception(f"The command {res.args} failed to run and returned the returncode={res.returncode}")
	res = subprocess.run(f"kmc_dump temp_union_max temp_dump_max; cat temp_dump_max | cut -f 2 | paste -sd+ - | bc",
						 shell=True, capture_output=True)
	if res.returncode != 0:
		raise Exception(f"The command {res.args} failed to run and returned the returncode={res.returncode}")
	wji_denominator = int(res.stdout)
	print(wji_denominator)

	# rm intermediate files
	os.remove("temp_intersect_min.kmc_pre")
	os.remove("temp_intersect_min.kmc_suf")
	os.remove("temp_union_max.kmc_pre")
	os.remove("temp_union_max.kmc_suf")
	if not keep_intermediate:
		os.remove("temp_dump_max")
		os.remove("temp_dump_min")

	# return WJI
	return wji_numerator * 1.0 / wji_denominator
# Ground truth object (independent of CEs)
class GroundTruth_WJI(object):
	"""
	Object for GT WJI calculation
	Ignore the read training_db because its for input file paths only
	Work for a whole range here
	"""
	def __init__(self, k_range: str, seq_file=None, rev_comp=True, temp_dir=None):
		self.rev_comp = rev_comp
		self.k_sizes = parsenumlist(k_range)
		self.seq_file = seq_file
		self.temp_dir = temp_dir
		if not os.path.exists(temp_dir):
			os.mkdir(temp_dir)

	# make a output file path by temp_dir + basename of input file
	def __kmc_output_name_converter(self, input_file: str, k_size: str) -> str:
		temp_dir = self.temp_dir
		return f"{os.path.join(temp_dir, os.path.basename(input_file))}_k_{k_size}"






####### other functions that may need to be added:
# insert into database: https://github.com/dkoslicki/CMash/blob/1485a0cde4685c52ae8cfb8bdaa5d1bad77eaac3/CMash/MinHash.py#L506
# delete from database: https://github.com/dkoslicki/CMash/blob/1485a0cde4685c52ae8cfb8bdaa5d1bad77eaac3/CMash/MinHash.py#L478
# merge 2 databases: https://github.com/dkoslicki/CMash/blob/1485a0cde4685c52ae8cfb8bdaa5d1bad77eaac3/CMash/MinHash.py#L563



######### Tests
def test_wji():
	E1 = CountEstimator(n=0, ksize=21)
	E2 = CountEstimator(n=0, ksize=21)

	E1._mins = [1,2,4,7]
	E2._mins = [1,2,5,6]

	E1._counts = [1,2,3,4]
	E2._counts = [1,1,2,2]

	assert E1.est_weighted_jaccard(E2) ==  2/8.0

def test_truncation():
	E3 = CountEstimator(n=0, ksize=4, rev_comp=False)
	E3._kmers = ['AAAA', 'AAAT', 'AAAC', 'AAGG']
	E3._counts = [1,2,3,4]
	# truncate to k=3
	E3.brute_force_truncation(3)
	assert E3._kmers == ['AAA', 'AAG']
	assert E3._counts == [6,4]
	# truncate to k=2
	E3.brute_force_truncation(2)
	assert E3._kmers == ['AA']
	assert E3._counts == [10]

def test_gt_wji():
	"""
	test the performance of ground truth WJI calculation
	:return:
	"""
	subprocess.run(f"echo TTAATTAA > test1.fq", shell=True)
	subprocess.run(f"echo AAATTTAA > test2.fq", shell=True)
	kmc_count("test1.fq", "out_test1", 4)
	kmc_count("test2.fq", "out_test2", 4)
	test_wji = kmc_intersection_count_for_wji("out_test1", "out_test2")
	os.remove("test1.fq")
	os.remove("test2.fq")
	print(test_wji)
	assert test_wji == 2 / 8.0
	
	
	












































