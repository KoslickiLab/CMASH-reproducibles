date

date
### local variables
[ -z $threads ] && threads=16
ltime="/usr/bin/time -av -o temp_runLog"

temp=$(which conda)
conda_path=$(echo ${temp%/*bin/conda})
if [ -f ${conda_path}/etc/profile.d/conda.sh ]; then
        . ${conda_path}/etc/profile.d/conda.sh
else
        echo "ERROR: conda path can't be corrected identified!!!"
        exit 1
fi
conda activate /data/sml6467/projects/202104_CMash_WABI_reproducibility/CMASH-reproducibles/src/conda_env/CMash_env_py37
python --version

### start running
file_name=random_30_genomes.txt
python /data/sml6467/github/CMASH-reproducibles/src/CMash_python_wrapper_of_CMash_github_to_perform_prefix_tree_bias_analysis.py -q ${file_name} -r ${file_name} -o True -t ${threads}  
cp est_JI_k60.csv trunc_JI_k60.csv


echo "Pipe done"
date

