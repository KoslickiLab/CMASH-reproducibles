# functions not used in the comparison pipe
# directly follow the usage
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


