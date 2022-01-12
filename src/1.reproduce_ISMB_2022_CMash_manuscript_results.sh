# To reproduce WABI 2021 CMash results
# Usage: bash <pipe>
threads=$1


date
### local variables
[ -z $threads ] && threads=16
ltime="/usr/bin/time -av -o temp_runLog"
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd ${pipe_path}/../1_ISMB_2022/
echo "Producing CMash results"
echo "This step may take a long time, it might be better to use nohup"
# the trunc_CI can be separate out from Bias_factor step to fasten the process
sleep 5 
echo "Pipe start!"



### Step1. download data
time_tag=$(date +"%Y%m%d_%H")
mkdir -p CMash_out_${time_tag}
cd CMash_out_${time_tag}
mkdir final_output
out_dir=$PWD
final=${out_dir}/final_output
mkdir -p ${final}/fig_1f_2_input
mkdir -p ${final}/fig_1f_3_input



date
echo "1. Downloading genome files!!!!!!"
sleep 2
# download 30 Brucella data
date
echo "1.1 Downloading 30 Brucella data"
for file in $(cut -f 3 ../data/Brucella_genus_30_genomes.txt); do
	suffix=$(echo ${file##*/})
	wget -q ${file}/${suffix}_genomic.fna.gz 2>/dev/null
done
mkdir Brucella_30 \
	&& mv GCA*_genomic.fna.gz ./Brucella_30 \
	&& cd ./Brucella_30 \
	&& readlink -f GCA*_genomic.fna.gz > ../fullpath_Brucella_30.txt \
	&& cd ..

date
# download 1000 random genomes
echo "1.2 Downloading 1000 random species"
for file in $(cut -f 3 ../data/random_1k_genomes.txt); do
	suffix=$(echo ${file##*/})
	wget -q ${file}/${suffix}_genomic.fna.gz 2>/dev/null
done
mkdir random_1000 \
	&& mv GCA*_genomic.fna.gz ./random_1000 \
	&& cd ./random_1000 \
	&& readlink -f GCA*_genomic.fna.gz > ../fullpath_random_1000.txt \
	&& cd ..

date
# download 200 genomes (from the 1000) for metagenomic simulation
echo "1.3 Downloading 200 genomes for simulation"
for file in $(cut -f 3 ../data/200_genomes_out_of_1000_for_simulation.txt); do
	suffix=$(echo ${file##*/})
	wget -q ${file}/${suffix}_genomic.fna.gz 2>/dev/null
done
mkdir simulation_200 \
	&& mv GCA*_genomic.fna.gz ./simulation_200 \
	&& cd ./simulation_200 \
	&& gunzip GCA*_genomic.fna.gz \
	&& cat GCA*_genomic.fna  > ../merged_200_genomes.fa \
	&& cd ..





### Step2. BBMap simulation and CI estimation for 1k genomes
# activate conda
temp=$(which conda)
conda_path=$(echo ${temp%/*bin/conda})
if [ -f ${conda_path}/etc/profile.d/conda.sh ]; then
        . ${conda_path}/etc/profile.d/conda.sh
else
        echo "ERROR: conda path can't be corrected identified!!!"
        exit 1
fi
unset conda_path
conda activate ${pipe_path}/conda_env/CMash_env_py37

# BBMap simulation
date
echo "2.1 BBMap simulation"
for depth in $(seq 2 2 20); do
	randomreads.sh ref=merged_200_genomes.fa out=BB_simulated_meta_${depth}m.fq reads=${depth}m len=150 metagenome
done

# CI estimation for all genomes
mkdir fig3_CI_estimation \
	&& mv BB_simulated_meta_*fq fig3_CI_estimation \
	&& cp fullpath_random_1000.txt fig3_CI_estimation \
	&& cd fig3_CI_estimation 

date
echo "2.2 Getting CI estimation for all 1000 ref genomes in metagenomic reads"
work_fig3=$PWD
ref=$(readlink -f fullpath_random_1000.txt)
temp_range=15-60-5
for file in $(ls BB_simulated_meta_*fq); do
	meta=$(readlink -f ${file})
	mkdir CMash_output_${file}
	cd CMash_output_${file}
	# CMash CI
	${ltime} MakeStreamingDNADatabase.py -k 60 -t ${threads} -n 2000 ${ref} TrainingDB_k60.h5
	mv temp_runLog build_ref.timelog
	${ltime} StreamingQueryDNADatabase.py -c 0 -l 0 -t ${threads}  --sensitive ${meta} TrainingDB_k60.h5 trunc_CI_result.csv ${temp_range}
	mv temp_runLog  get_CI.timelog
	mkdir trunc_CI_out
	mv TrainingDB_k60* trunc_CI_out
	mv *timelog trunc_CI_out
	mv trunc_CI_result.csv trunc_CI_out
	# Est CI
	for i in $(seq 15 5 60); do
		## build ref for each k
		${ltime} MakeStreamingDNADatabase.py -k ${i} -t ${threads} -n 2000 ${ref} TrainingDB_k${i}.h5
		mv temp_runLog build_ref_k${i}.timelog
		${ltime} StreamingQueryDNADatabase.py -c 0 -l 0 -t ${threads}  --sensitive ${meta} TrainingDB_k${i}.h5 est_CI_k${i}_result.csv ${i}-${i}-1
		mv temp_runLog  get_CI_k${i}.timelog
	done
	mkdir est_CI_out
	mv *timelog est_CI_out
	mv TrainingDB_k* est_CI_out
	mv est_CI_k*result.csv est_CI_out
	cd ${work_fig3}
done

date
echo "2.3 Collect time/space usage for CI estimation"
echo -e "depth\tspace\ttotal_time\tref_time\trunning_time" > time_summary_trunc_CI.txt
echo -e "depth\tspace\ttotal_time\tref_time\trunning_time" > time_summary_est_CI.txt

get_usr_sys_time () {
	prefix=$1 #use prefix*suffix to accept multipe inputs
	suffix=$2
	# mutltiple files will contain extra ":", so use the ")" as first layer filter and then use ":"
	sys=$(grep "System" ${prefix}*${suffix} | cut -d")" -f 2 | cut -d":" -f 2 | cut -d" " -f 2 | awk '{s+=$1}END{print s}')
	usr=$(grep "User" ${prefix}*${suffix} | cut -d")" -f 2 | cut -d":" -f 2 | cut -d" " -f 2 | awk '{s+=$1}END{print s}')
	sum=$(echo ${sys} + ${usr} | bc )
	echo ${sum}
}

for i in $(seq 2 2 20); do
	out_depth="${i}m"
	out_name_suffix="CMash_output_BB_simulated_meta_${i}m.fq"
	cd ${out_name_suffix}
	# infor from trunc_CI
	cd trunc_CI_out
	trunc_ci_ref_time=$(get_usr_sys_time build_ref timelog) #sec
	trunc_ci_proc_time=$(get_usr_sys_time get_CI timelog)
	trunc_ci_ref_size=$(ls -al *tst | awk '{s+=$5}END{print s}') #byte
	trunc_ci_ref_size=$(echo ${trunc_ci_ref_size} / 1024 / 1024 | bc ) #MB
	trunc_ci_total_time=$(echo ${trunc_ci_ref_time} + ${trunc_ci_proc_time} | bc )
	# infor from est_CI
	cd ../est_CI_out
	est_ci_ref_size=$(ls -al *tst | awk '{s+=$5}END{print s}') #byte
	est_ci_ref_size=$(echo ${est_ci_ref_size} / 1024 / 1024 | bc ) #MB
	est_ci_ref_time=$(get_usr_sys_time build_ref timelog) #sec
	est_ci_proc_time=$(get_usr_sys_time get_CI timelog)
	est_ci_total_time=$(echo ${est_ci_ref_time} + ${est_ci_proc_time} | bc )
	# summary result
	cd ${work_fig3}
	echo -e "${out_depth}\t${trunc_ci_ref_size}\t${trunc_ci_total_time}\t${trunc_ci_ref_time}\t${trunc_ci_proc_time}" >> time_summary_trunc_CI.txt
	echo -e "${out_depth}\t${est_ci_ref_size}\t${est_ci_total_time}\t${est_ci_ref_time}\t${est_ci_proc_time}" >> time_summary_est_CI.txt
done

for i in $(seq 2 2 20); do
	out_depth="${i}m"
	out_name_suffix="CMash_output_BB_simulated_meta_${i}m.fq"
	cd ${out_name_suffix}
	cp trunc_CI_out/trunc_CI_result.csv ${work_fig3}/trunc_CI_results_${out_depth}.csv
	cd est_CI_out/
	### need to merge est_CI output (pre-sorted)
	for k in $(seq 15 5 60); do
		cp est_CI_k${k}_result.csv  ${work_fig3}/est_CI_results_${out_depth}_k${k}.csv
	done
	cd ${work_fig3}
done
mkdir CI_outputs
mv trunc_CI_results_*csv ./CI_outputs
mv est_CI_results_*csv ./CI_outputs
cp ./CI_outputs/*csv  ${final}/fig_1f_3_input
cp time_summary_est_CI.txt ${final}/fig_1f_3_input
cp time_summary_trunc_CI.txt ${final}/fig_1f_3_input
# for fig3a, the cumusum of ref size
cd $(ls -d CMash_output_* | head -1) \
	&& ls -al trunc_CI_out/*tst | awk '{print $5}' > ../cum_space_trunc_CI.txt \
	&& ls -al est_CI_out/*tst | awk '{print $5}' > ../cum_space_est_CI.txt \
	&& cd .. \
	&& cp cum_space_*txt ${final}/fig_1f_3_input




### Step3, pairwise JI from containment MinHash for 30 genomes in Brucella genus
date
echo "3. get pairwise JI for 30 genomes in Brucella genus"
cd ${out_dir}
mkdir fig2_JI_estimation
cp fullpath_Brucella_30.txt ./fig2_JI_estimation
cd fig2_JI_estimation
file_name=$(readlink -f fullpath_Brucella_30.txt)
python  ${pipe_path}/CMash_python_wrapper_of_CMash_github_to_perform_prefix_tree_bias_analysis.py -q ${file_name} -r ${file_name} -o True -t ${threads}
find . -name "*JI*csv" | xargs cp -t ${final}/fig_1f_2_input
# trunc_JI_k60 is equivalent to est_JI_k60 because no truncation heppened
cd ${final}/fig_1f_2_input \
	&& cp est_JI_k60.csv trunc_JI_k60.csv






### Step4, genrating figures from output files
date
echo "4. generate figures"
cd ${final}
python ${pipe_path}/RECOMB_2021_generate_plots.py








date
echo "Whole pipe done!"




