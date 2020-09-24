### Implement WTST trial in bash
# Usage: bash <script> 

### Pipe infor:
# v1.00, last update 09/24/2020
#
###



### set variables
while getopts q:r:k:c:t:n:h opts
do
	case "$opts" in
		q) query="$OPTARG";;		# query file
		r) ref="$OPTARG";;		# ref file
		k) maxk="$OPTARG";;		# max size K of k-mer
		c) input_range="$OPTARG";;	# range of size to check, format: start-end-gap
		t) threads="$OPTARG";;		# thread number for CMash
		n) num_hashes="$OPTARG";;	# hash function number, default 2000
		h) echo "
Run WTST trial to compare estimates from truncation and brute force calculation.

Usage: bash <pipe> -q <query> -r <ref> -k <max_k> -c <range> 

Options:
-q: a txt file containing the path of input query genome files (fa or fa.gz)
-r: a txt file containing the path of input ref genome files (fa or fa.gz)
-k: max size of k-mer to use
-c: range of the truncation, e.g. 10-60-5
-t: number of threads to use, default 8
-n: number of hashes to use in MinHash estimates, default 2000
"
exit;;
[?]) echo "use -h for help"
exit;;
esac
done

date
echo "Pipe start......"
echo "Checking variables......"

# check required parameters
if [ -z "$query" ] || [ -z "$ref" ] || [ -z "$maxk" ] || [ -z "$input_range" ]; then
	echo "Missing required parameters!!!"
	exit 1
fi

# set default
[ -z "$threads" ] && threads=8
[ -z "$num_hashes" ] && num_hashes=2000

# range adjustment
# a range should start from at least 2
# if a range can't be divided by the gap, the start value would be adjusted, and the max value would be kept
temp_range=`echo $input_range | awk -F"-" '{ if(($1==1)) print $1+1"-"$2"-"$3; else print $1"-"$2"-"$3}'`
  r_start=`echo $temp_range | cut -d"-" -f 1`
  r_end=`echo $temp_range | cut -d"-" -f 2`
  r_gap=`echo $temp_range | cut -d"-" -f 3`
  r_adj_start=$((r_start+(r_end-r_start)%r_gap))
temp_range=${r_adj_start}-${r_end}-${r_gap}

# local variables
ltime="/usr/bin/time -av -o temp_runLog"
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
echo "Input checking success!!!"
echo " "

# active conda
echo "Activate conda env......"
temp=$(which conda)
conda_path=$(echo ${temp%/*bin/conda})
if [ -f ${conda_path}/etc/profile.d/conda.sh ]; then
	. ${conda_path}/etc/profile.d/conda.sh
else
	echo "Conda path can't be corrected identified......"
	exit 1
fi
conda activate ${pipe_path}/py37_WTST && echo "Conda activated!!!" || (echo "Failed to activate conda env......" && exit 1)

# transform path to absolute path
query=$(readlink -f $query)
ref=$(readlink -f $ref)



### build k-mer sketches, BF, and WTST for input files
${ltime} python ${pipe_path}/create_ref.py -q ${query} -r ${ref} -k ${maxk} -n ${num_hashes} -t ${threads}




### read database to do comparison
${ltime} python ${pipe_path}/compute_similarity.py -q <query_WTST>  -r <ref_WTST>  -c ${range} -t ${threads}



### brute-forcely compute the CI (no truncation)
for i in <range>; do
	${ltime} python ${pipe_path}/create_ref.py -q ${query} -r ${ref} -k ${i} -n ${num_hashes} -t ${threads}
	${ltime} python ${pipe_path}/compute_similarity.py -q <query_WTST_i> -r <ref_WTST_i> -c <i-i-1> -t ${threads}
done

bash ${pipe_path}/brute_force_result_collection.sh



### compare the results
python ${pipe_path}/compare_results.py <trunc_result> <brute_force_result>



date
echo "pipe done"
















