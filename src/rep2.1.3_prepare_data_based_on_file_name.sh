# for auto simulation:
# the meta file must be fasta (unzip)
# the ref can be just gz
# all be absolute path
# run in the current folder of the file

meta=$1
ref=$2

if [ -z $meta ] || [ -z $ref ]; then
	echo "missing input"
fi

workdir=$PWD

# clean file
cat $meta | sed '/^$/d' | sed 's/estimated_CI_k61_//g' | sed 's/"//g' >  cleaned_meta.txt
cat $ref | sed '/^$/d' | sed 's/"//g' > cleaned_ref.txt

# grep file
mkdir fasta_meta
find /data/sml6467/resource/NCBI_Genbank_bacteria_genome_20200525/NCBI_Genbank_bacteria_genome/genome_files | grep -w -f cleaned_meta.txt | xargs cp -t fasta_meta
cd fasta_meta
gunzip *
find $PWD -name "*fna" > ../path_meta.txt
cd ..

find /data/sml6467/resource/NCBI_Genbank_bacteria_genome_20200525/NCBI_Genbank_bacteria_genome/genome_files | grep -w -f cleaned_ref.txt > path_ref.txt

echo "bingo"
date


