# CMash-reproducibles
For reproducibility analysis of [CMash](https://github.com/dkoslicki/CMash) related tasks:
- [1. WABI 2021 CMash manuscript](#1_wabi_2021)
- [2. CMash full paper](#2_cmash_full)


## Install dependencies
1. Install conda or [miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) (version >= 4.6)
2. Clone this repo
```
git clone https://github.com/KoslickiLab/CMASH-reproducibles.git
```
3. Get required tools/resources (only need to do once)
```
cd CMASH-reproducibles/src
bash 0.install_all_required_dependency_run_once.sh
```
4. Uninstall everything if necessary
```
cd CMASH-reproducibles/src
bash uninstall.sh
```



## 1. WABI 2021 CMash manuscript <a name="1_wabi_2021"></a>
This is to reproduce the results in the [CMash manuscript for WABI 2021](https://www.overleaf.com/project/6074afe80ace05a51a7e71e9).
Please follow **Install dependencies** above first to install all required dependencies.
1. regenerate all the output data (may take more than 1 day)
```
cd CMASH-reproducibles/src
nohup bash 1.reproduce_WABI_2021_CMash_manuscript_results.sh  &  #accept 1 positional parameter for thread number (default 16)
```
2. find the results
```
cd CMASH-reproducibles/1_WABI_2021_CMash_manuscript
ls -d CMash_out_*  #output folder: CMash_out_${time_tag}
```
3. folder structure
```
# final_output: stores all final output files and Fig 1f, 2, and 3
# fig2_JI_estimation: intermediate outputs for pairwise JI estimation within Brucella genus
# fig3_CI_estimation: intermediate outputs for containment estimation of 1000 random genomes in the simulated metagenomic data
# Brucella_30 / random_1000 / simulation_200: downloaded genome files
```





