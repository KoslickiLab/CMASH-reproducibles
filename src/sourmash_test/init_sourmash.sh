# Setup sourmash env to run analysis + timing
# Manually run!!!
# last update: 08/17/2020


### create Env 
conda create -y -n Sourmash_env_py37 python=3.7
conda activate Sourmash_env_py37
conda install -y -c conda-forge -c bioconda sourmash
conda install -y -c bioconda khmer

### test sourmash
sourmash compute -h
conda deactivate


