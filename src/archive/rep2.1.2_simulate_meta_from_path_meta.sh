# pre-requsite:
# 1. the meta file must be unzipped
# 2. need the taxid-species file 
# 3. a demo mini_confir.ini for CAMISIM

taxid_species_file="/data/sml6467/projects/202002_CMash_test/results/rep2_CI_change_against_gradient_ref/simulation/2_runCAMISIM/taxid-species_id-link.txt"   # obtained from NCBI note file
demo_config="/data/sml6467/projects/202002_CMash_test/results/rep2_CI_change_against_gradient_ref/simulation/2_runCAMISIM/demo_mini_config.ini"


out_name=$1 # the out name
input_file=$2 #the full path to "path_meta.txt" file
ref_file=$3
if [ -z $out_name ] || [ -z $input_file ] || [ -z $ref_file ]; then
	echo "please sepcify name, input_file, and ref file"
	exit 1
fi



cd /data/sml6467/projects/202002_CMash_test/results/rep2_CI_change_against_gradient_ref/simulation/auto_run_simulation
mkdir output_${out_name}
cd output_${out_name}

awk '{print "genome"NR"\t"$1}' ${input_file} > genome_to_id.tsv
echo -e "genome_ID\tOTU\tNCBI_ID\tnovelty_category" > metadata.tsv

while read p;
do
aa=$(echo $p | awk '{print $1}')
temp=$(echo $p | awk '{print $2}')
bb=$(wc -l $temp | awk '{print $1}')
name=$(echo ${temp##*/})
name=$(echo ${name%_genomic.fna})
temp2=$(grep $name $taxid_species_file)
cc=$(echo $temp2 | awk '{print $2}')
echo $temp2 >> record.txt
echo -e "$aa\t$bb\t$cc\tKnown_strain" >> metadata.tsv
done <  genome_to_id.tsv

# run CAMISIM
genome_num=$(wc -l $input_file | awk '{print $1}')
cat $demo_config | sed "s/replace_genome_total/$genome_num/g" | sed "s/replace_genome_real/$genome_num/g" > mini_config_${out_name}.ini
cp genome_to_id.tsv /data/sml6467/software/CAMISIM/shaopeng_testrun
cp metadata.tsv /data/sml6467/software/CAMISIM/shaopeng_testrun
cp mini_config_${out_name}.ini /data/sml6467/software/CAMISIM/shaopeng_testrun
workdir=$PWD
cd /data/sml6467/software/CAMISIM
conda_path="/data/sml6467/software/miniconda3"
. ${conda_path}/etc/profile.d/conda.sh
conda activate py27_CAMISIM
python metagenomesimulation.py shaopeng_testrun/mini_config_${out_name}.ini
conda deactivate
# mv results back
mv 


cd $workdir
# run BBMap
for file in $(cat $input_file); do
	echo $file
	cat $file >> merged_all.fa
done

# current 10M >1 coverage, so didn't use that parameter
bash /data/sml6467/software/bbmap/randomreads.sh ref=merged_all.fa out=BBMap_simulated_meta_10M.fq reads=10m len=150 metagenome 

# run autotest pipe on the simulated data and ref
readlink -f BBMap_simulated_meta_10M.fq > bb_path.txt
query_bbmap=$(readlink -f bb_path.txt)

mkdir run_comparison
cd run_comparison
conda activate CMash-env

bash /data/sml6467/projects/202002_CMash_test/src/rep1.1_meta_again_ref_need_polish.sh -q ${query_bbmap}   -r ${ref_file} 


