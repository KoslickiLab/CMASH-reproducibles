#!/bin/bash

###### Pipeline information:
# v1.0, last update 07/08/2020
# 1 or more metagenomes vs ref file for GroundTruth / Estimated / Truncated CIs

######################################################################################################
###### Setup variables
while getopts q:r:k:c:g:d:t:n:e:h opts
do
	case "$opts" in
		q) query="$OPTARG";;		# query file
		r) ref="$OPTARG";;			# ref file
		k) maxk="$OPTARG";;			# max size of k-mer
		c) input_range="$OPTARG";;		# range of size to check, format: start-end-gap
		g) CMash="$OPTARG";;			# path to the CMash repo
		d) conda_path="$OPTARG";;	# conda path
		t) threads="$OPTARG";;		# thread number
		n) num_hashes="$OPTARG";;	# hash function number
		e) cmash="$OPTARG";;		# cmash environment to use
		h) echo "
Calculate the containment index (CI) change between query data and ref data with specified k-mer size.

Usage: bash <pipe.sh> -q <query> -r <ref> -k <max_k> -c <kmer_range> -g <CMash_repo> -d <conda_path>

Options:
-q: query file, a txt file with the absolute path of the input metagenome data
-r: ref file, a txt file absolute path of all the reference genome data
-k: max size of k-mer used in the analysis
-c: range of k-mer size to check, must in start-end-gap format, the starting poing would be adjusted if the gap is not a divisor of the range,
	e.g. 10-30-5 means sizes of k-mer to check are 10,15,20,25,30
-g: absolute path to the github repo "CMASH"
-d: absolute path to the conda folder, this is used to activate conda env inside bash script
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
if [ -z "$query" ] || [ -z "$ref" ] || [ -z "$maxk" ] || [ -z "$input_range" ] || [ -z "$CMash" ] || [ -z "$conda_path" ] || [ -z "$cmash" ]; then
echo "Missing input parameter!!!"
	exit 1
fi
query=$(readlink -f $query)
ref=$(readlink -f $ref)
conda_path=$(readlink -f $conda_path)
CMash=$(readlink -f $CMash)

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
. ${conda_path}/etc/profile.d/conda.sh  
conda activate $cmash
ltime="/usr/bin/time -av -o temp_runLog"
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
summary_code="${pipe_path}/rep1_summary.py"
# test CMash python env match
export PYTHONPATH=${CMash}:$PYTHONPATH
testFile=$(python -c "from CMash import MinHash as MH; print(MH.__file__)")
correctFile=${CMash}/CMash/MinHash.py
[ "$testFile" == "$correctFile" ] && echo "Files are correct" || exit 1
scriptDir=${CMash}/scripts
moduleDir=${CMash}/CMash
echo "Checking done, start to process file"




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
	${ltime} python ${scriptDir}/MakeStreamingDNADatabase.py -k ${i} -v -t ${threads} -n ${num_hashes}  ${ref} TrainingDB_k${i}.h5
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
	${ltime} python ${scriptDir}/StreamingQueryDNADatabase.py ${file} ../ref_db_${temp_range}/TrainingDB_k${maxk}.h5 truncation_${name}_results.csv ${temp_range} -v -c 0 -l 0 -t ${threads}  --sensitive
	mv temp_runLog truncation_${name}_results.log
done
cd ..

# Step3, estimated CI and ground truth CI
for i in $(seq ${r_adj_start} ${r_gap} ${maxk})
do
	echo "running on k ${i}"
	for file in `cat ${query}`
	do
		echo $file
		name=`echo ${file##*/}`
		# estimated CI
		${ltime} python ${scriptDir}/StreamingQueryDNADatabase.py ${file} ./ref_db_${temp_range}/TrainingDB_k${i}.h5 estimated_CI_k${i}_${name}_results.csv ${i}-${i}-1 -v -c 0 -l 0 -t ${threads} --sensitive 
		mv temp_runLog estimated_CI_k${i}_${name}_results.log
		# ground truth CI
		${ltime} python ${moduleDir}/GroundTruth.py ${file} ./ref_db_${temp_range}/TrainingDB_k${i}.h5  ground_truth_CI_k${i}_${name}_results.csv ${i}-${i}-1 -c 0 -l 0 
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

conda deactivate
echo "whole pipe done"
date


