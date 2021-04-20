#!/bin/bash

###### Pipeline information:
# v1.01, last update 07/19/2020
# 1 or more metagenomes vs ref file for GroundTruth / Estimated / Truncated CIs
# now need to activate the CMash Env in advance to run this pipe!!!
# 
# Update history:
# V1.01: put all possible dependencies into conda env


######################################################################################################
###### Setup variables
while getopts q:r:k:c:t:n:h opts
do
	case "$opts" in
		q) query="$OPTARG";;		# query file
		r) ref="$OPTARG";;			# ref file
		k) maxk="$OPTARG";;			# max size of k-mer
		c) input_range="$OPTARG";;		# range of size to check, format: start-end-gap
		t) threads="$OPTARG";;		# thread number
		n) num_hashes="$OPTARG";;	# hash function number
		h) echo "
Calculate the containment index (CI) change between query data and ref data with specified k-mer size.

Usage: bash <pipe.sh> -q <query> -r <ref> -k <max_k> -c <kmer_range> -g <CMash_repo> -d <conda_path>

Options:
-q: query file, a txt file with the absolute path of the input metagenome data
-r: ref file, a txt file absolute path of all the reference genome data
-k: max size of k-mer used in the analysis
-c: range of k-mer size to check, must in start-end-gap format, the starting poing would be adjusted if the gap is not a divisor of the range,
	e.g. 10-30-5 means sizes of k-mer to check are 10,15,20,25,30
-t: number of threads to use for CMash analysis, default 48
-n: number of hashes functions to use in CMash, default 2000
-h: help information
"
exit;;
[?]) echo "use -h for help"
exit;;
esac
done

### check input parameter
if [ -z "$query" ] || [ -z "$ref" ] || [ -z "$maxk" ] || [ -z "$input_range" ]; then
	echo "Missing input parameter!!!"
	exit 1
fi
query=$(readlink -f $query)
ref=$(readlink -f $ref)
[ -z "$threads" ] && threads=48
[ -z "$num_hashes" ] && num_hashes=2000

### range adjustment
temp_range=`echo $input_range | awk -F"-" '{ if(($1==1)) print $1+1"-"$2"-"$3; else print $1"-"$2"-"$3}'`
  r_start=`echo $temp_range | cut -d"-" -f 1`
  r_end=`echo $temp_range | cut -d"-" -f 2`
  r_gap=`echo $temp_range | cut -d"-" -f 3`
  r_adj_start=$((r_start+(r_end-r_start)%r_gap))
temp_range=${r_adj_start}-${r_end}-${r_gap}

### dependent parameters
ltime="/usr/bin/time -av -o temp_runLog"
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
summary_code="${pipe_path}/rep1_summary.py"



###### Running CMash
time_tag=`date +"%m_%d_%H-%M"`
mkdir CMash_output_${time_tag}
cd CMash_output_${time_tag}

# Step1, build the ref database for all ks in range
mkdir ref_db_${temp_range}
cd ref_db_${temp_range}
for i in $(seq ${r_adj_start} ${r_gap} ${maxk})
do
	echo "generating reference for k=${i}"
	${ltime} MakeStreamingDNADatabase.py -k ${i} -v -t ${threads} -n ${num_hashes}  ${ref} TrainingDB_k${i}.h5
        mv temp_runLog TrainingDB_k${i}.log	2> /dev/null
done
cd ..

# Step2, truncated CI for all k-point
mkdir truncated_CI
cd truncated_CI
for file in `cat ${query}`
do
	echo $file
	name=`echo ${file##*/}`
	${ltime} StreamingQueryDNADatabase.py ${file} ../ref_db_${temp_range}/TrainingDB_k${maxk}.h5 truncation_${name}_results.csv ${temp_range} -v -c 0 -l 0 -t ${threads}  --sensitive
	mv temp_runLog truncation_${name}_results.log
done
cd ..

# Step3, estimated CI and ground truth CI
### currently GT.py hasn't been added to the PATH
temp=$(which StreamingQueryDNADatabase.py)
temp=$(readlink -f $temp)
env_dir=$(echo ${temp%/bin/StreamingQueryDNADatabase.py})
gt_py=$(find $env_dir -name "GroundTruth.py" | grep "CMASH_Env_py37" | head -1) # incease multiple matches during code adjustment

for i in $(seq ${r_adj_start} ${r_gap} ${maxk})
do
	echo "running on k ${i}"
	for file in `cat ${query}`
	do
		echo $file
		name=`echo ${file##*/}`
		# estimated CI
		${ltime} StreamingQueryDNADatabase.py ${file} ./ref_db_${temp_range}/TrainingDB_k${i}.h5 estimated_CI_k${i}_${name}_results.csv ${i}-${i}-1 -v -c 0 -l 0 -t ${threads} --sensitive 
		mv temp_runLog estimated_CI_k${i}_${name}_results.log
		# ground truth CI
		${ltime} python ${gt_py} ${file} ./ref_db_${temp_range}/TrainingDB_k${i}.h5  ground_truth_CI_k${i}_${name}_results.csv ${i}-${i}-1 -c 0 -l 0 
		mv temp_runLog ground_truth_CI_k${i}_${name}_results.log
	done
done

mkdir ground_truth_CI
mv ground_truth_CI_k* ground_truth_CI
mkdir estimated_CI
mv estimated_CI_k* estimated_CI


# Step4, summaize results by py
mkdir summary
python ${summary_code} ${temp_range} ${query}

echo "whole pipe done"
date


