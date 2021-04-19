### install necessary dependencies for WTST test
# only need to run once

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



### build Conda env 
conda create -y -p ${pipe_path}/py37_WTST python=3.7
conda activate ${pipe_path}/py37_WTST
conda install -y -c bioconda cmash
conda install -y -c anaconda seaborn
conda install -y -c bioconda khmer
conda install -y -c bioconda kmc
conda deactivate


echo "Pipe done"
date


