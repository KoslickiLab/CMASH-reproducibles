# ISMB 2022 CMash input data

There are three files for genomes used in the manuscript:
```
1. "Brucella_genus_30_genomes.txt" contains the 30 genomes used for pairwise JI in Fig2
2. "random_1k_genomes.txt" and 
3. "200_genomes_out_of_1000_for_simulation.txt" are the 1000 random genomes used for CI estimation in Fig3
```
The format of the 3 files are: taxid, sci name, download link.


In addition, 4 groups (of genomes) from the 1000 random genomes were picked to explore the distribution
 of bias factor:
```
1. "sourmash_1kgenome_k30.csv": the pairwise JI values bewteen these 1000 genomes
2. "rand1k_group1_JI_0.3-0.6.txt": name of a small group with moderate pairwise JI values
3. "rand1k_group1_JI_0.5-8.txt": name of a small group with moderate-high pairwise JI values
4. "rand1k_group1_JI_0.9-1.txt": name of a small group with high pairwise JI values
5. "rand1k_group1_JI_0.txt": a random sample of 10 genomes with validated 0 (or near 0) JI values
```




