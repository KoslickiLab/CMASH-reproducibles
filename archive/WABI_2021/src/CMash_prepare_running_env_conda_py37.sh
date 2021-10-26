# prepare conda env (py3.7) for CMash analysis (note: py3.8 is not compatible due to pkg conflict)
# the conda env will be installed in "src/conda_env"
# usage: bash <pipe>



# to call conda inside a bash
temp=$(which conda)
conda_path=$(echo ${temp%/*bin/conda})
if [ -f ${conda_path}/etc/profile.d/conda.sh ]; then
        . ${conda_path}/etc/profile.d/conda.sh
else
        echo "ERROR: conda path can't be corrected identified!!!"
        exit 1
fi
unset conda_path



# location
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
conda_path=${pipe_path}/conda_env
git_path=${pipe_path}/git_repo



# install CMash env
cd ${conda_path}
conda create -y -p ${PWD}/CMash_env_py37 python=3.7
conda activate ${PWD}/CMash_env_py37
conda install -y -c bioconda cmash
conda install -y -c anaconda seaborn
conda install -y -c bioconda kmc
conda install -y -c bioconda taxonkit
conda install -y -c bioconda insilicoseq
conda install -y -c bioconda bbmap
conda deactivate



# download NCBI taxdump for taxonkit
cd CMash_env_py37
mkdir -p taxonkit_ncbi
cd taxonkit_ncbi
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz && rm taxdump.tar.gz



# download CMash git
cd ${git_path}
git clone https://github.com/dkoslicki/CMash.git



date
echo "pipe done"
