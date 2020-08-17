# Running cmd for sourmash
# Last update 08/17/2020

date
### dir
cd /data/sml6467/projects/202002_CMash_test/results/try_sourmash
conda activate Sourmash_env_py37
sourmash compute -h 
ltime="/usr/bin/time -av -o temp_runLog"

### Demo code from the tool paper 
# Large-scale sequence comparisons with sourmash
# Use local files for testing
cp /data/sml6467/projects/202002_CMash_test/results/20200727_test_V102/CMASH-reproducibles/task1_k-mer_truncation/example/output_07_27_16-50/raw_query/GCA_001022035.1_ASM102203v1_genomic.fna . 
cp /data/sml6467/projects/202002_CMash_test/results/20200727_test_V102/CMASH-reproducibles/task1_k-mer_truncation/example/output_07_27_16-50/raw_query/GCA_011765425.1_ASM1176542v1_genomic.fna .

query="GCA_001022035.1_ASM102203v1_genomic.fna"
ref="GCA_011765425.1_ASM1176542v1_genomic.fna"

# 1.1 hash signature (no k-mer infor, only hash)
# also works for fastq file
# may need trim-low-abund.py to trim low freq reads for metagenome
${ltime} sourmash compute -k 21,31,51 \
	--scaled 2000 \
	--track-abundance \
	-o query.sig  ${query}
mv temp_runLog query_sig.log


sourmash compute -k 21,31,51 \
        --scaled 2000 \
        --track-abundance \
        -o ref.sig  ${ref}
### rename to track the changes if necessary


# 1.2 compare sigs
${ltime} sourmash compare -k 31 --csv test.csv *.sig 
mv temp_runLog compare_sig.log
### Jaccard for unweighted k-mer; Cosine for weighted k-mer


# 1.3 MDS clustering analysis by R

# 1.4 Tetranucleotide Frequency Clustering
### logic: build signature by k=4, then call pairwise comparison
mkdir tetra
cd tetra
wget https://osf.io/xusfa/download -O Oe6.scaffolds.sub.fa.gz

sourmash compute -k 4 \
	--scaled 1 \
	--track-abundance \
	--singleton \
	-o Oe6.scaffolds.sub.comp \
	Oe6.scaffolds.sub.fa.gz

### ssems a numpy issue need to be fixed:
# https://github.com/automl/auto-sklearn/issues/667
sourmash plot --labels \
	--vmin .4 \
	Oe6.scaffolds.sub.comp

cd ..

# 1.5 Best-first search (SBTMH all contaiment)
mkdir search
cd search
### create ecoli database from signatures
mkdir escherichia–sigs 
cd escherichia–sigs
wget https://osf.io/pc76j/download -O escherichia–sigs.tar.gz
tar xzf escherichia–sigs.tar.gz
rm escherichia–sigs.tar.gz
cd ..
sourmash index -k 31 ecolidb ./escherichia–sigs/ecoli-*sig

### get ecoli reads
wget https://osf.io/26xm9/download -O ecoli–reads.khmer.fq.gz

sourmash compute -k 31 \
	--scaled 2000 \
	-o ecoli-reads.sig \
	ecoli–reads.khmer.fq.gz

### search
sourmash search -k 31 \
	ecoli-reads.sig ecolidb \
	--containment

# 1.6 Best-first gather (iteratively SBTMH substract best hit) 
sourmash gather -k 31 \
	ecoli-reads.sig ecolidb





