## Notes

### For "1_WABI_2021_..." folder:
That's incomplete version for method part (no error analysis).

This is to reproduce the results in the [CMash manuscript for RECOMB 2021](https://www.overleaf.com/project/61666320f8392d53f75c1135).
Please follow **Install dependencies** above first to install all required dependencies.

(Put the folder under main repo rolder)  
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
 



