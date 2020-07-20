#!/bin/bash

###### Pipeline information:
# v1.01, last update 07/19/2020. 
# Install dependencies for CMash reproducible analysis
# Only need to run once
######



date
### conda path: put this cmd out to switch to standard example 
temp=$(which conda) 
conda_path=$(echo ${temp%/*bin/conda})
if [ -f ${conda_path}/etc/profile.d/conda.sh ]; then
	. ${conda_path}/etc/profile.d/conda.sh
else
	echo "ERROR: conda path can't be corrected identified!!!"
	exit 1
fi



### get paths for pipe
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
repo=$(echo ${pipe_path%/src})



### TEMP: add the py38 CMash to run "GroundTruth.py"
conda env create -p ${pipe_path}/temp_CMash_py38 -f ${pipe_path}/CMash_env.yml
conda activate ${pipe_path}/temp_CMash_py38
conda install -y -c bioconda cmash
conda deactivate


### build CMash conda env to run CMash analysis
conda create -y -p ${pipe_path}/CMASH_Env_py37 python=3.7
conda activate ${pipe_path}/CMASH_Env_py37
conda install -y -c bioconda cmash
conda install -y -c anaconda seaborn
conda deactivate



### build CAMISIM (must be py27) to run CAMISIM + BBMAP
conda create -y -p ${pipe_path}/Simulation_Env_py27 python=2.7
conda activate ${pipe_path}/Simulation_Env_py27
conda install -y -c bioconda samtools=1.9
conda install -y -c bioconda bbmap
conda install -y -c anaconda biopython
conda deactivate



### download the NCBI GenBank bacteria database
wget -O ${repo}/src/NCBI_GenBank_bacteria_assembly_summary.txt  ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
cd ${repo}/src
awk -F '\t'  '{if($12=="Complete Genome" && $11=="latest") print $20}' NCBI_GenBank_bacteria_assembly_summary.txt > NCBI_GenBank_download_link.txt



### download CAMISIM
git clone https://github.com/CAMI-challenge/CAMISIM.git
# add new taxdamp file 
### this may not always work due to the structure change of the original file
### but CAMISIM doesn't take the original input, so I've to update it
wget -O CAMISIM/tools/ncbi-taxonomy_20200705.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
cd CAMISIM/tools/
mkdir NCBI
mv ncbi-taxonomy_20200705.tar.gz ./NCBI
cd NCBI
tar -xzf ncbi-taxonomy_20200705.tar.gz && rm ncbi-taxonomy_20200705.tar.gz
cd ..
tar -zcvf ncbi-taxonomy_20200705.tar.gz NCBI



echo "Pipe done"
date
