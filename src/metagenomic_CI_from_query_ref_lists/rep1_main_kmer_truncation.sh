#!/bin/bash

###### Pipeline information:
# v1.02, last update 07/29/2020
# Taking 2 files containing the absolute path for geneome data as query / ref data, estimate the CI change between GroundTruth CI, Estimated CI and Truncated CI
# Input parameters:
#	1. query file: names of species from NCBI GenBank database, e.g. GCA_002786755.1_ASM278675v1
#	2. ref file: names of species from NCBI GenBank database, e.g. GCA_002786755.1_ASM278675v1
#	3. maxk: the maximum size of k
#	4. range: series of k-mer length to check the CI change
#	5. repo: absolute path to the github repo
#	6. conda path: absolute path to the conda main folder (activate env inside bash)
# It will create a new folder named "output_M_D_H-M" in the current dir
#
#
# Update history:
# V1.01: put all possible dependencies into conda env; 
#		 add uninstall.sh;
# V1.02: add depth parameter for BBMap (coverage) and CAMISIM (output size / GB)




######################################################################################################
###### Setup variables
date 
while getopts q:r:k:c:s:d:t:h opts
do
	case "$opts" in
		q) query="$OPTARG";;		# query file
		r) ref="$OPTARG";;			# ref file
		k) maxk="$OPTARG";;			# max size of k-mer
		c) range="$OPTARG";;		# range of size to check, format: start-end-gap
		s) size="$OPTARG";;			# depth for CAMISIM (output file size / GB), BBMap coverage, default 2
		d) conda_path="$OPTARG";;	# path to conda
		t) threads="$OPTARG";;		# thread number for CMash
		h) echo "
Calculate the containment index (CI) change between query data and ref data with specified k-mer size.

Usage: bash <pipe.sh> -q <query_file> -r <ref_file> -k <max_k> -c <kmer_range>

Options:
-q: query file, a txt file with each line containing the name of one COMPLETE genome file from GenBank, 
	e.g. GCA_002786755.1_ASM278675v1
-r: ref file, a txt file with each line containing the name of one COMPLETE genome file from GenBank, 
	e.g. GCA_002786755.1_ASM278675v1
-k: max size of k-mer used in the analysis
-c: range of k-mer size to check, must in start-end-gap format, the starting poing would be adjusted if the gap is not a divisor of the range,
	e.g. 10-30-5 means sizes of k-mer to check are 10,15,20,25,30
-d: path to the main conda folder (not the bin folder), e.g. ~/user/miniconda3 
-t: number of threads to use for CMash analysis, default 48
-h: help information
"
exit;;
[?]) echo "use -h for help"
exit;;
esac
done



pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
repo=$(echo ${pipe_path%/src})



### temporary usage to check the main conda path
if [ -z "$conda_path" ]; then
	temp=$(which conda) 
	conda_path=$(echo ${temp%/*bin/conda})
	if [ -f ${conda_path}/etc/profile.d/conda.sh ]; then
		. ${conda_path}/etc/profile.d/conda.sh
	else
		echo "ERROR: conda path can't be corrected identified!!!"
		exit 1
	fi
fi

### check depth
[ -z "$size" ] && size=2


### check input parameter
if [ -z "$query" ] || [ -z "$ref" ] || [ -z "$maxk" ] || [ -z "$range" ] || [ -z "$repo" ]; then
	echo "Missing input parameter!!!"
	exit 1
fi
[ -z "$threads"] && threads=48


### convert path to absolute path
query=$(readlink -f $query)
ref=$(readlink -f $ref)
conda_path=$(readlink -f $conda_path)


### read prepared information
gb_ref=${repo}/src/NCBI_GenBank_bacteria_assembly_summary.txt
gb_download=${repo}/src/NCBI_GenBank_download_link.txt


### other var
ltime="/usr/bin/time -av -o temp_runLog"



######################################################################################################
###### Data preparation
# output dir
time_tag=`date +"%m_%d_%H-%M"`
mkdir output_${time_tag}
cd output_${time_tag}
date > running_record.log
echo -e "\nInput files:\nQuery\t$query\nRef\t$ref\nGithub_repo\t$repo\nConda_path\t$conda_path\nmaxk\t$maxk\nrange\t$range\ncoverage(depth)\t$size\n\n" >> running_record.log
workdir=${PWD}


# download query genomes, while CAMISIM can't handle gz file, unzip all fna.gz
mkdir raw_query
cd raw_query
cat $query | sed '/^$/d' | sed "s/['\"]//g" | sed 's/_genomic.fna.gz//g'  > cleaned_query.txt
for file in $(cat cleaned_query.txt)
do
  echo $file
  wget $(grep $file ${gb_download})/${file}_genomic.fna.gz
done
gunzip *_genomic.fna.gz 	
find ${PWD} -name "*_genomic.fna" > ../query_path.txt
cd ..


# download ref genomes
mkdir raw_ref
cd raw_ref
cat $ref | sed '/^$/d' | sed "s/['\"]//g" | sed 's/_genomic.fna.gz//g'  > cleaned_ref.txt
for file in $(cat cleaned_ref.txt)
do
  echo $file
  wget $(grep $file ${gb_download})/${file}_genomic.fna.gz
done
find ${PWD} -name "*_genomic.fna.gz" > ../ref_path.txt
cd ..


### Run simulation based on BBMAP / CAMISIM
conda activate ${pipe_path}/Simulation_Env_py27
# BBMap simulation
mkdir BBMap_simu
for file in $(cat query_path.txt); do
  cat $file >> BBMap_simu/merged_all.fa 
done
cd BBMap_simu
${ltime} randomreads.sh ref=merged_all.fa out=BBMap_simulated_meta_3x.fq coverage=${size} len=150 metagenome  &>  running_record_BBMap.log
readlink -f BBMap_simulated_meta_3x.fq > ../bb_meta_path.txt
cd ..


# CAMISIM simulation
mkdir CAMISIM_simu
awk '{print "genome"NR"\t"$1}' query_path.txt > genome_to_id.tsv
genome_num=$(wc -l query_path.txt | awk '{print $1}')
mv genome_to_id.tsv CAMISIM_simu
cd CAMISIM_simu

### prepare input files: metadata.tsc and genome_to_id.tsv
echo -e "genome_ID\tOTU\tNCBI_ID\tnovelty_category" > metadata.tsv
while read p;
do
	aa=$(echo $p | awk '{print $1}')
	temp=$(echo $p | awk '{print $2}')
	bb=$(wc -l $temp | awk '{print $1}')
	name=$(echo ${temp##*/})
	name=$(echo ${name%_genomic.fna*})
	cc=$(grep $name ${gb_ref} | awk -F '\t'  '{if($12=="Complete Genome" && $11=="latest") print $7}')
	grep $name ${gb_ref} | awk -F '\t'  '{if($12=="Complete Genome" && $11=="latest") print $0}' >> record.txt
	echo -e "$aa\t$bb\t$cc\tKnown_strain" >> metadata.tsv
done <  genome_to_id.tsv

### prepare input file: the config file
cp ${repo}/src/CAMISIM/defaults/mini_config.ini new_miniconfig.ini
mv genome_to_id.tsv new_genome_to_id.tsv
mv metadata.tsv new_metadata.tsv
# replace corresponding values in the config file
cami_size=$(awk "BEGIN {print ${size}/10}")
sed -i "s/size=0.1/size=$cami_size/g" new_miniconfig.ini
sed -i 's/metadata.tsv/new_metadata.tsv/g' new_miniconfig.ini
sed -i 's/genome_to_id.tsv/new_genome_to_id.tsv/g' new_miniconfig.ini
sed -i "s/genomes_total=24/genomes_total=$genome_num/g" new_miniconfig.ini
sed -i "s/genomes_real=24/genomes_real=$genome_num/g" new_miniconfig.ini
sed -i 's/ncbi-taxonomy_20170222.tar.gz/ncbi-taxonomy_20200705.tar.gz/g' new_miniconfig.ini
sed -i 's#tools/samtools-1.3/samtools#../Simulation_Env_py27/bin/samtools#g' new_miniconfig.ini
cp new* ${repo}/src/CAMISIM/defaults/
cd ${repo}/src/CAMISIM/
${ltime} python metagenomesimulation.py defaults/new_miniconfig.ini  &> running_record_CAMISIM.log 
mv temp_runLog ${workdir}/CAMISIM_simu
mv running_record_CAMISIM.log ${workdir}/CAMISIM_simu
conda deactivate
if [ -d "out" ]; then
	mv out ${workdir}/CAMISIM_simu
	cd ${workdir}/CAMISIM_simu/out
	for file in $(find . -name "anonymous_reads.fq.gz"); do
		cat $file >> merged_CAMISIM_anonymous.fq.gz
	done
	readlink -f merged_CAMISIM_anonymous.fq.gz > ../../cami_meta_path.txt
fi


###### Run CMash
# while the running task for metagenome against ref is an independent task and would be used frequently
# it is an independent bash script and called here
cd ${workdir}
conda activate ${pipe_path}/CMASH_Env_py37
bash ${repo}/src/rep1_meta_vs_ref_CI_compare.sh -q bb_meta_path.txt -r ref_path.txt -k ${maxk} -c ${range}  -t ${threads}  &> CMash_BBMap.log
mv CMash_output_* Output_CMash_BBMap
bash ${repo}/src/rep1_meta_vs_ref_CI_compare.sh -q cami_meta_path.txt -r ref_path.txt -k ${maxk} -c ${range} -t ${threads} &> CMash_CAMISIM.log
mv CMash_output_* Output_CMash_CAMISIM
conda deactivate

date
date >> running_record.log
echo "whole pipe done"

  
