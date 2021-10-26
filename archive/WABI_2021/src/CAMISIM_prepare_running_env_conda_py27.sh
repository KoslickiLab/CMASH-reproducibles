# prepare conda env (py2.7) for usage of CAMISIM
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

# install cami
cd ${git_path}
git clone https://github.com/CAMI-challenge/CAMISIM
### update taxdamp file
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir NCBI
mv taxdump.tar.gz ./NCBI/
cd NCBI
tar -xzf taxdump.tar.gz && rm taxdump.tar.gz
cd ..
time_tag=$(date +"%Y%m%d")
tar -zcvf ncbi-taxonomy_${time_tag}.tar.gz NCBI && rm -r NCBI
mv ncbi-taxonomy_${time_tag}.tar.gz ./CAMISIM/tools/
### create a sample "mini_config.ini" file with updated taxdump
cd ./CAMISIM/defaults
mkdir personal_run && cd personal_run
cp ../mini_config.ini templete_mini_config.ini
sed -i "s/ncbi-taxonomy_20170222.tar.gz/ncbi-taxonomy_${time_tag}.tar.gz/g" templete_mini_config.ini
sed -i 's/max_strains_per_otu=1/max_strains_per_otu=100/g' templete_mini_config.ini
sed -i 's#tools/samtools-1.3/samtools#../../conda_env/cami_env_py27/bin/samtools#g' templete_mini_config.ini



# prepare conda env for CAMISIM
cd ${conda_path}
conda create -y -p ${PWD}/cami_env_py27 python=2.7
conda activate ${PWD}/cami_env_py27
conda install -y -c bioconda samtools=1.9
conda install -y -c bioconda bbmap
conda install -y -c anaconda biopython
conda deactivate



date
echo "Finish building CAMI conda env"
