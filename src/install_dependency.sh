#!/bin/bash

###### Pipeline information:
# Install dependencies for CMash reproducible analysis
# Only need to run once
######



date
###### read parameters
while getopts p:h opts
do
	case "$opts" in
		p) conda_env_path="$OPTARG";;	# optional, conda env storage location if needed
		h) echo "
Install dependencies for CMASH-reproducible analysis.

Usage: bash <pipe.sh>

Options:
		-p: optional, specify the location to store the Conda environment
		-h: help information
"
exit;;
[?]) echo "use -h for help"
exit;;
esac
done

pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
repo=$(echo ${pipe_path%/src})


###### Build Conda Env
# in case the channel is missing
conda config --add channels conda-forge

# write env path 
date > ${repo}/src/source.txt 

# build env
if [ -z "$conda_env_path" ]; then
	# Install CMash-Env
	echo "Building environment for CMash"
	conda env create -f ${repo}/src/CMash_env.yml \
	&& echo -e "Env_for_CMash\tCMash_env" >> ${repo}/src/source.txt
	# Install CAMISIM-Env
	echo "Building environment for CAMISIM"
	conda env create -f ${repo}/src/CAMISIM_env.yml \
	&& echo -e "Env_for_CAMISIM\tCAMISIM_env" >> ${repo}/src/source.txt
else
	# Install CMash-Env
	echo "Building environment for CMash"
	conda env create -p ${conda_env_path}/CMash-env -f ${repo}/src/CMash_env.yml \
	&& echo -e "Env_for_CMash\t${conda_env_path}/CMash-env" >> ${repo}/src/source.txt
	# Install CAMISIM-Env
	echo "Building environment for CAMISIM"
	conda env create -p ${conda_env_path}/CAMISIM-env -f ${repo}/src/CAMISIM_env.yml \
	&& echo -e "Env_for_CAMISIM\t${conda_env_path}/CAMISIM-env" >> ${repo}/src/source.txt
fi
echo " " >> ${repo}/src/source.txt



###### Download the NCBI GenBank bacteria database
wget -O ${repo}/src/NCBI_GenBank_bacteria_assembly_summary.txt  ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
cd ${repo}/src
awk -F '\t'  '{if($12=="Complete Genome" && $11=="latest") print $20}' NCBI_GenBank_bacteria_assembly_summary.txt > NCBI_GenBank_download_link.txt


###### Download CMash
git clone https://github.com/dkoslicki/CMash.git


###### Prepare CAMISIM and BBMap
# CAMISIM
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

# BBMAP 38.86
cd ${repo}/src
wget —referer https://sourceforge.net/projects/bbmap/files/BBMap_38.86.tar.gz/download 
mv download BBMap_38.86.tar.gz \
	&& tar -xvzf  BBMap_38.86.tar.gz \
	&& rm BBMap_38.86.tar.gz \
	&& rm index.html 2>/dev/null
	



echo "Pipe done"
date
