import time
import os
import screed
import re
import bisect
import khmer
import h5py
import numpy as np
import multiprocessing
from multiprocessing import Pool


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
    :param mins1:
    :param mins2:
    :param counts1:
    :param counts2:
    :return:
    """
    sum_min = 0 #for numerator
    sum_max = 0 #for denominator
    i = 0
    j = 0
    processed = 0
    # if a MH value appears only in one list, then min(W(A), W(B)) = 0
    # so we only count the other one if there is no match
    try: #avoid the outofbound error in last round
        while processed <= min(len(mins1), len(mins2)):
            while mins1[i] < mins2[j]:
                sum_max += counts1[i] #W(B)=0, so add W(A)
                i += 1
                processed += 1
            while mins1[i] > mins2[j]:
                sum_max += counts2[j] #W(A)=0, so add W(B)
                j += 1
                processed += 1
            if mins1[i] == mins2[j]:
                #skip the value = max_prime check becaue their weights are all 0, so min/max are 0s
                sum_min += min(counts1[i], counts2[j])
                sum_max += max(counts1[i], counts2[j])
                i += 1
                j += 1
    except IndexError:
        WJI = sum_min / float(sum_max)
        return WJI

# calculate standard JI from previous WJI function by truncting the counts==0
def not_necessary_ji(mins1, mins2, counts1, counts2):
    i = len(counts1)-1
    j = len(counts2)-1
    while counts1[i] == 0:
        mins1.pop()
        i -= 1
    while counts2[j] == 0:
        mins2.pop()
        j -= 1

    i = 0
    j = 0
    processed = 0
    sum_min = 0 #for numerator
    sum_max = 0 #for denominator
    # repeat the function above
    try: #avoid the outofbound error in last round
        while processed <= min(len(mins1), len(mins2)):
            while mins1[i] < mins2[j]:
                sum_max += 1 #W(B)=0, so add W(A)
                i += 1
                processed += 1
            while mins1[i] > mins2[j]:
                sum_max += 1 #W(A)=0, so add W(B)
                j += 1
                processed += 1
            if mins1[i] == mins2[j]:
                #skip the value = max_prime check becaue their weights are all 0, so min/max are 0s
                sum_min += 1
                sum_max += 1
                i += 1
                j += 1
    except IndexError:
        JI = sum_min / float(sum_max)
        return JI

# parallel wrapper for WJI calculation
def unwrap_jaccard_vector(arg):
    return arg[0].est_weighted_jaccard(arg[1])

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
    def __init__(self, n=None, max_prime=None, ksize=None, input_file_name=None, hash_list=None, rev_comp=True):
        if n is None:
            raise Exception
        if max_prime is None:
            raise Exception
        if ksize is None:
            raise Exception
        if input_file_name is None:
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
                    self.add(kmer, rev_comp)
        else:
            for sub_seq in seq_split_onlyACTG:
                if len(sub_seq)>=self.ksize:        #in case of small chunk
                    self.add_sequence(sub_seq, rev_comp)

    # add a kmer into the sketch
    def add(self, kmer, rev_comp):
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
        if h >= _mins[-1]:
            return

        # insert kmer into the sketch
        i = bisect.bisect_left(_mins, h)  # find index to insert h
        if _mins[i] == h: #already in sketch
            _counts[i] += 1
        else:
            #h not in sketch, insert
            _mins.insert(i, h)
            _mins.pop()
            _counts.insert(i, 1)
            _counts.pop()
            _kmers.insert(i, kmer)
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
























