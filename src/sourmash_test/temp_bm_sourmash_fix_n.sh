# Running cmd for sourmash
# Last update 08/23/2020

### Setup variable
date
while getopts q:r:c:n:t:d:ah opts
do
	case "$opts" in
		q) query="$OPTARG";;		# query file containing all abs path
		r) ref="$OPTARG";;		# ref file containing all abs path
		c) range="$OPTARG";;		# range of size to check, format: start-end-gap
		n) num_hash="$OPTARG";;		# the num_hash to use, default 2000
		t) threads="$OPTARG";;		# number of threads to use, default 48
		a) abundance="yes";;		# enable abundance tracking, default is no
		d) conda_path="$OPTARG";;	# conda path to use
		h) echo "
Benchmarking for sourmash
Usage: bash <script> -q <query> -r <ref> -c <range>
"
exit;;
[?]) echo "use -h for help"
exit;;
esac
done
################################



### var validation
# check input
if [ -z "$query" ] || [ -z "$ref" ] || [ -z "$range" ] || [ -z "$conda_path" ]; then
	echo "Missing input parameter!!!"
	exit 1
fi
[ -z "$num_hash" ] && num_hash=2000
[ -z "$threads" ] && threads=48
[ -z "$abundance" ] && abun_indicator="" || abun_indicator="--track-abundance" # for abundance option in building signature

ltime="/usr/bin/time -av -o temp_runLog"
query=$(readlink -f $query)
ref=$(readlink -f $ref)
conda_path=$(readlink -f $conda_path)
# range adjustment
temp_range=`echo $range | awk -F"-" '{ if(($1==1)) print $1+1"-"$2"-"$3; else print $1"-"$2"-"$3}'`
  r_start=`echo $temp_range | cut -d"-" -f 1`
  r_end=`echo $temp_range | cut -d"-" -f 2`
  r_gap=`echo $temp_range | cut -d"-" -f 3`
  r_adj_start=$((r_start+(r_end-r_start)%r_gap))
temp_range=${r_adj_start}-${r_end}-${r_gap}
k_mer_sets=$(seq -s ","  $r_adj_start $r_gap $r_end) # for sourmash usage


### local var, can be deleted if the env is merged into the CMash-env
conda_path="/data/sml6467/software/miniconda3/"
. ${conda_path}/etc/profile.d/conda.sh
conda activate Sourmash_env_py37
#cd /data/sml6467/projects/202002_CMash_test/results/try_sourmash
time_tag=`date +"%m_%d_%H-%M"`
mkdir output_${time_tag}
cd output_${time_tag}
###


### Build sigunature of all ref files
### rename to track the changes if necessary
# the threads option works for 10x bam files only
${ltime} sourmash compute -k ${k_mer_sets} \
	-n ${num_hash} -p ${threads} ${abun_indicator} \
	-o ref.sig $(cat $ref | paste -s)
mv temp_runLog record_sourmash_compute_ref.log

### Build sig for all query files 1 by 1 (because each input is one study)
# if put all sigs together, there are unnecessary pairwise estimate between 2 query files
# may need trim-low-abund.py to trim low freq reads for metagenome
for file in $(cat $query); do
	echo $file
	name=`echo ${file##*/}`
	${ltime} sourmash compute -k ${k_mer_sets} \
		-n ${num_hash} -p ${threads} ${abun_indicator} \
		-o ${name}.sig ${file}
	mv temp_runLog record_sourmash_compute_${name}.log
	# do the similarity est at the meantime
	for k in $(seq $r_adj_start $r_gap $r_end); do
		${ltime} sourmash compare -k ${k} -p  ${threads}  --csv output_${name}_k${k}.csv ${name}.sig ref.sig
		mv temp_runLog record_sourmash_compare_${name}_k${k}.log
	done
done



### grab all the time
mkdir time_summary
cd time_summary
echo "collecting all times"
date > time_summary.txt
echo -e "Time record (User + Sys) for sourmash cmds  \n" >> time_summary.txt
echo -e "Signature creationg time\t" >> time_summary.txt


# signature creation time
for file in $(ls ../record_sourmash_compute_*); do
	echo $file
	name=$(echo ${file#../record_sourmash_compute_})
	name=$(echo ${name%.log})
	echo $name
	temp_time=$(grep "User" ${file} | head -1 | cut -d":" -f 2 | cut -d" " -f 2)
	echo -e "User_${name}\t${temp_time}" >> part1_sig_creationg_time.txt
	temp_time=$(grep "System" ${file} | head -1 | cut -d":" -f 2 | cut -d" " -f 2)
	echo -e "Sys_${name}\t${temp_time}" >> part1_sig_creationg_time.txt
done

# JI est time
for file in $(ls ../record_sourmash_compare_*); do
	echo $file
	name=$(echo ${file#../record_sourmash_compare_})
	name=$(echo ${name%.log})
	echo $name
	temp_time=$(grep "User" ${file} | head -1 | cut -d":" -f 2 | cut -d" " -f 2)
	echo -e "User_${name}\t${temp_time}" >> part2_similarity_est_time.txt
	temp_time=$(grep "System" ${file} | head -1 | cut -d":" -f 2 | cut -d" " -f 2)
	echo -e "Sys_${name}\t${temp_time}" >> part2_similarity_est_time.txt
done

# write into summary
### ref sig time
t_sig_ref=$(grep -e "User_ref" -e "Sys_ref" part1_sig_creationg_time.txt | awk '{s+=$2}END{print s}')
echo -e "All ref files\t${t_sig_ref}" >> time_summary.txt
### each query sig time
for file in $(cat $query); do
	name=`echo ${file##*/}`
	t_sig_q=$(grep -e "User_${name}" -e "Sys_${name}" part1_sig_creationg_time.txt | awk '{s+=$2}END{print s}')
	echo -e "${name}\t${t_sig_q}" >> time_summary.txt
done
### each query compare time
echo -e "\nSimilarity estimation time\t" >> time_summary.txt
for file in $(cat $query); do
        name=`echo ${file##*/}`
        t_sig_q=$(grep -e "User_${name}" -e "Sys_${name}" part2_similarity_est_time.txt | awk '{s+=$2}END{print s}')
	echo -e "${name}\t${t_sig_q}" >> time_summary.txt
done


echo "whole pipe done"
conda deactivate
date












