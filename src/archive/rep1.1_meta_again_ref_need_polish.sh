#!/bin/bash
date
echo "pipe start"

### read parameters
while getopts r:q:o:c:h opts
do case "$opts" in
	r) ref_file="$OPTARG";; # ref database path
	q) query_file="$OPTARG";; # query path
	o) outdir="$OPTARG";; # output dir
	c) CMash="$OPTARG";; # the CMash locol copy to use
	h) echo "
Description:
!!! Havan't fully polished !!!
Compare the ground truth / estimated / truncated containment index of arbitrary input metagenome and references.

Usage: bash <pipe> -r <ref> -q <query>

Optional parameters:
-c <cmash_dir>: locol copy of CMash to use, default the master branch
"
exit;;
[?]) echo "Use -h for help information"
exit;;
esac
done


### check parameters
if [ -z $ref_file ] || [ -z $query_file ]; then
	echo "Please specify both ref and query"
	exit 1
fi

if [ -z "$outdir" ]; then
        outdir=$PWD
	echo "outdir is $outdir"
fi

if [ -z "$CMash" ]; then
       CMash="/data/sml6467/github/CMash_master"
fi

### input parameters to be set
input_range="10-60-5"
conda_path="/data/sml6467/software/miniconda3"
python_exec="python3.8"
threads=48
num_hashes=2000

### parameter check
# check input range: 1) can't start from 1; 2) maxk should be covered
[ -z "$input_range" ] && echo "Please specify your testing range in the format start-end-gap" && exit 1
temp_range=`echo $input_range | awk -F"-" '{ if(($1==1)) print $1+1"-"$2"-"$3; else print $1"-"$2"-"$3}'`
  r_start=`echo $temp_range | cut -d"-" -f 1`
  r_end=`echo $temp_range | cut -d"-" -f 2`
  r_gap=`echo $temp_range | cut -d"-" -f 3`
  r_adj_start=$((r_start+(r_end-r_start)%r_gap))
temp_range=${r_adj_start}-${r_end}-${r_gap}
maxk=${r_end}


### dependent paramebers (no need to touch)
# activate conda env inside bash
. ${conda_path}/etc/profile.d/conda.sh  
conda activate CMash-env
# simplify time cmd
ltime="/usr/bin/time -av -o temp_runLog"
# dir of pipe folder
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
summary_code="${pipe_path}/rep1_summary.py"
# test CMash python env match
export PYTHONPATH=${CMash}:$PYTHONPATH
testFile=$(python -c "from CMash import MinHash as MH; print(MH.__file__)")
correctFile=${CMash}/CMash/MinHash.py
[ "$testFile" == "$correctFile" ] && echo "Files are correct" || exit 1
scriptDir=${CMash}/scripts
moduleDir=${CMash}/CMash
time_tag=`date +"%m_%d_%H-%M"`
echo "Checking done, start to process file"

### Running the py code
cd ${outdir} 
mkdir output_${time_tag}
cd output_${time_tag}

# Step1, build the ref database for all ks in range
mkdir ref_db_${temp_range}
cd ref_db_${temp_range}
for i in $(seq ${r_adj_start} ${r_gap} ${maxk})
do
	echo "generating reference for k=${i}"
	${ltime} ${python_exec} ${scriptDir}/MakeStreamingDNADatabase.py -k ${i} -v -t ${threads} -n ${num_hashes}  ${ref_file} TrainingDB_k${i}.h5
        mv temp_runLog TrainingDB_k${i}.log	2> /dev/null
done
cd ..


# Step2, truncated CI for all k-point
### temp adjust to test code
mkdir truncated_CI
cd truncated_CI
for file in `cat ${query_file}`
do
	echo $file
	name=`echo ${file##*/}`
	${ltime} ${python_exec} ${scriptDir}/StreamingQueryDNADatabase.py ${file} ../ref_db_${temp_range}/TrainingDB_k${maxk}.h5 truncation_${name}_results.csv ${temp_range} -v -c 0 -l 0 -t ${threads}  --sensitive
	mv temp_runLog truncation_${name}_results.log
done
cd ..

# Step3, estimated CI and ground truth CI
for i in $(seq ${r_adj_start} ${r_gap} ${maxk})
do
	echo "running on k ${i}"
	for file in `cat ${query_file}`
	do
		echo $file
		name=`echo ${file##*/}`
		# estimated CI
		${ltime} ${python_exec} ${scriptDir}/StreamingQueryDNADatabase.py ${file} ./ref_db_${temp_range}/TrainingDB_k${i}.h5 estimated_CI_k${i}_${name}_results.csv ${i}-${i}-1 -v -c 0 -l 0 -t ${threads} --sensitive 
		mv temp_runLog estimated_CI_k${i}_${name}_results.log
		# ground truth CI
		${ltime} ${python_exec} ${moduleDir}/GroundTruth.py ${file} ./ref_db_${temp_range}/TrainingDB_k${i}.h5  ground_truth_CI_k${i}_${name}_results.csv ${i}-${i}-1 -c 0 -l 0 
		mv temp_runLog ground_truth_CI_k${i}_${name}_results.log
	done
done

mkdir ground_truth_CI
mv ground_truth_CI_k* ground_truth_CI
mkdir estimated_CI
mv estimated_CI_k* estimated_CI


# Step3, summaize results by py
mkdir summary
${python_exec} ${summary_code} ${temp_range} ${query_file}

conda deactivate
date
