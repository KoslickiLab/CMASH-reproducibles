# CMash-reproducibles
For reproducibility of [Cmash](https://github.com/dkoslicki/CMash) tasks/issues:
1. [K-mer truncation](https://github.com/dkoslicki/CMash/issues/20)

## Install dependencies  
1. [miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
2. CMash environment  
Git clone CMash +  [yml file](https://github.com/dkoslicki/CMash/blob/master/tests/CMash-env-Ubuntu-18.05.yml) + seaborn
3. CAMISIM environment (for non-Docker usage)
Git clone CAMISIM + py27 conda env + biopython + biom-format

## Reproducibility
1. To reproduce reuslts at [xxx](), please:  
   - [go to folder]() and download the data by `bash xxx`
   - run the cmd `bash xxxxx`

## Task1: how would truncation of long k-mers affect containment index?
Comparing the containment index generated by truncated k-mers (prefix of original long k-mers) and original k-mers of a large k (such as 60).

#### Input file:
1. ref file: a file containing absolute path of all reference data in analysis
2. query file: a file containing absolute path of all query metagenome data in analysis (must be unzipped, as CAMISIM can't work with zipped file)

#### Output file:
1. CAMISIM and BBMap simulated reads in fastq format
2. Tenery search tree of a range of k-mer database
3. Estimated / ground_truch / truncated containment index (stored in csv file) of the input metagenome files
4. Summary plot of the comparison

#### Example
[Click here](https://github.com/KoslickiLab/CMASH-reproducibles/tree/master/task1_K-mer_truncation/example) for an example input and output.  
Input:
1. ref file is "path_ref.txt"
2. query file is "path_meta.txt"

Output:
1. estimated_CI: store the estimated containment index between metagenome and reference datasets
2. ground_truth_CI: store the ground truch containment index between metagenome and reference datasets
3. truncated_CI: store the containment index estimated by truncated k-mer
4. ref_df: store the tenery search tree of each k size (intentially left empty here as the file size is large)
5. summary: the summary plot



