import os
import pandas as pd
import numpy as np
import seaborn as sb
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import math
import statistics



####################### define functions
# read data
### clean matrix for symmatric input: pairwise JI matrix is symmatric so we only need half of the matrix
def upper_triangle(input_df, lower_tril=True):
	if lower_tril:
		print("Ploting by the lower triangle of the input matrix!!!")
		out_df = input_df.where(np.tril(np.ones(input_df.shape), -1).astype(np.bool))
	else:
		out_df = input_df.where(np.triu(np.ones(input_df.shape), 1).astype(np.bool))
	return out_df



### For fig2 collect data and store in a dict: gt_CI, trunc_CI
### use a dict to store all csv for gt_CI / trunc_CI
def read_all_k_df(sect_name, k_range, symmatric=True, maxk_for_bf=None):
	"""
	read dfs for all k values of each section and store in a dict
	:param sect_name: est, gt, trunc, bias
	:param k_range: a list of k values, e.g. [15, 20, 25]
	:param symmatric: if the matrix is symmatric, only lower half will be kept
	:param maxk_for_bf: maxk is needed for bias factor usage
	:return: a dict storing all dfs
	"""
	
	if sect_name == "est":
		prefix = "est_JI_k"
	elif sect_name == "gt":
		prefix = "GroundTruth_JI_k"
	elif sect_name == "trunc":
		prefix = "trunc_JI_k"
	elif sect_name == "bias":
		prefix = "bias_factor_k"
		if maxk_for_bf == None:
			raise Exception("Please specify maxk when reading BF dfs")
	else:
		raise Exception("Sect_name must be one of ['est', 'gt', 'trunc', 'bias']")
	
	out_dict = dict()
	
	for k in k_range:
		if sect_name == "bias":
			f_name = prefix + str(k) + "_to_k" + str(maxk_for_bf) + ".csv"
		else:
			f_name = prefix + str(k) + ".csv"
		# read df for a single k value
		temp_df = pd.read_csv(f_name, header=0, index_col=0)
		if symmatric:
			temp_df = upper_triangle(temp_df)
		out_dict[str(k)] = temp_df
	
	return out_dict
### transfer dict to a df for plot purpose (will lose structure infor)
def dict_to_df_by_k(input_dict, k_range):
	# merge the values into a big DF
	p_out = pd.DataFrame()
	for k in k_range:
		temp_colname = 'k=' + str(k)
		all_values_1k = input_dict[str(k)].stack().to_numpy()
		p_out[temp_colname] = all_values_1k
	return p_out



### For fig3 collect data and return the consistency matrix for fig3a
def collect_ci_dif_data(k_range, name_key, cutoff=0.05):
	out_df = pd.DataFrame()
	for key in name_key:
		f_name = "trunc_CI_results_" + key + ".csv"
		trunc_df = pd.read_csv(f_name, header=0, index_col=0)
		temp_list = []
		for k in k_range:
			f_name = "est_CI_results_" + key + "_k" + str(k) + ".csv"
			temp_k_df = pd.read_csv(f_name, header=0, index_col=0)
			temp_trunc_df = trunc_df[["k=" + str(k)]]
			# remember to merge by colname: 2 files have different sort!!!
			temp_check = pd.concat([temp_k_df, temp_trunc_df], axis=1, join="inner")
			temp_list.append(sum(abs(temp_check.iloc[:, 0] - temp_check.iloc[:, 1]) <= cutoff) / 1000)
		out_df[str(key)] = temp_list
		out_df.index = k_range
	
	return out_df



### generate figures
### get dif matrix and plot fig1f.left
### Fig1f.left: pick ONLY 1 pair of data for a line plot
def get_dif_matrix_and_fig1f(d_trunc, d_gt, k_range):
	# merge out values into a big DF
	abs_out = pd.DataFrame()
	rela_out = pd.DataFrame()
	for k in k_range:
		temp_colname = 'k=' + str(k)
		# abs error
		dif_df = d_trunc[str(k)] - d_gt[str(k)]
		all_abs_values = dif_df.stack().to_numpy()
		abs_out[temp_colname] = all_abs_values
		# rela error
		rela_df = dif_df / (d_gt[str(k)] + 0.00001)
		all_rela_values = rela_df.stack().to_numpy()
		rela_out[temp_colname] = all_rela_values
		
	# generate fig1f.left
	all_keys = list(d_trunc.keys())
	line_tgt = (0, 1)  # just pick a pair for the single lineplot, can be arbitrary
	if math.isnan(d_trunc[all_keys[0]].iloc[line_tgt]):
		line_tgt = (1, 0)
	### collect 2 lists for the line plot
	line_trunc = []
	line_gt = []
	for ele in all_keys:
		line_trunc.append(d_trunc[ele].iloc[line_tgt])
		line_gt.append(d_gt[ele].iloc[line_tgt])
	### generate line plot in fig1f.left
	fig, axs = plt.subplots(1, 1)
	plt.rcParams.update({'font.size': 10})
	axs.set_ylim([0, 1])
	axs.set(ylabel="JI value", xlabel="k size")
	sb.lineplot(x=all_keys, y=line_trunc, ax=axs, color='#FFA500', marker="o", linewidth=3, markersize=15)
	sb.lineplot(x=all_keys, y=line_gt, ax=axs, color="black", marker="o", linewidth=3, markersize=15)
	legend_element = [
		Line2D([0], [0], marker='o', color='#FFA500', label='CMash JI', markerfacecolor='#FFA500', markersize=12),
		Line2D([0], [0], marker='o', color='black', label='GT JI', markerfacecolor='black', markersize=12)]
	axs.legend(handles=legend_element, loc='upper center', fontsize=22)
	fig.savefig("Fig1f_left_single_lineplot_Trunc_GT.png", dpi=200)
	plt.close(fig)
	
	# return the df for fig2
	return [abs_out, rela_out]


### Fig1f.right
def generate_fig1f_right(ref_size_list, trunc_time, est_time):
	trunc_size = ref_size_list[-1]
	est_size = sum(ref_size_list)
	trunc_ref_time = math.ceil(statistics.mean(list(trunc_time['ref_time']))/60)  # avg minute
	est_ref_time = math.ceil(statistics.mean(list(est_time['ref_time']))/60)  # avg sec
	
	# tune input vector
	sect = ["Truncation", "Classic MinHash"]
	ref_size = [trunc_size, est_size]
	ref_time = [trunc_ref_time, est_ref_time]
	
	# make barplot
	fig, axs = plt.subplots(1, 2, figsize=(12, 8))
	plt.rcParams.update({'font.size': 10})
	sb.barplot(x=sect, y=ref_size, palette=['#FFA500', 'black'], ax=axs[0])
	# axs[0].set_ylabel("MB", loc='top')
	axs[0].set_ylim([0, 1500])
	axs[0].set(xticklabels=[])
	axs[0].legend(title='Database size', title_fontsize=30, loc='upper center')
	for p in axs[0].patches:
		axs[0].annotate(str(round(p.get_height())) + " MB", (p.get_x() + p.get_width() / 2., p.get_height()), size=25,
		                ha='center', va='center', xytext=(0, 10), textcoords='offset points')
	sb.barplot(x=sect, y=ref_time, palette=['#FFA500', 'black'], ax=axs[1])
	axs[1].set_ylim([0, 2500])
	axs[1].set(xticklabels=[])
	axs[1].legend(title='CPU time', title_fontsize=30, loc='upper center')
	for p in axs[1].patches:
		axs[1].annotate(str(round(p.get_height())) + " min", (p.get_x() + p.get_width() / 2., p.get_height()), size=25,
		                ha='center', va='center', xytext=(0, 10), textcoords='offset points')
	fig.savefig("Fig1f_right_space_time_simple_compare.png", dpi=200)
	plt.close(fig)


### Fig2：
def generate_fig2(df_trunc, df_gt, df_abs_dif, df_rela_dif):
	# generate Fig2 from Trunc/GT processed matrix
	fig, axs = plt.subplots(2, 2, figsize=(20, 20))
	plt.rcParams.update({'font.size': 35})
	# ab: swarmplot of GT and Trunc JI
	sb.swarmplot(data=df_gt, ax=axs[0,0], color="black", size=1, linewidth=0.3)
	axs[0,0].set(ylabel="Ground truth Jaccard indices")
	axs[0,0].set_ylim([0, 1])
	sb.swarmplot(data=df_trunc, ax=axs[0, 1], color="black", size=1, linewidth=0.3)
	axs[0, 1].set(ylabel="CMash Jaccard indices")
	axs[0, 1].set_ylim([0, 1])
	# cd: absolute error and relative error
	sb.boxplot(data=df_abs_dif, ax=axs[1,0], color="black")
	axs[1,0].set(ylabel="CMash - Ground truth")
	axs[1,0].set_ylim([-0.3, 0.3])
	sb.boxplot(data=df_rela_dif, ax=axs[1,1], color="black")
	axs[1,1].set(ylabel="Relative error: CMash to Ground truth")
	axs[1,1].set_ylim([-0.3, 0.3])
	for i in [0,1]:
		for j in [0,1]:
			plt.setp(axs[i,j].get_xticklabels(), rotation=45)
	#fig.suptitle(plot_name)
	fig.savefig("Fig2_compare_trunc_gt_ji.png", dpi=200)
	plt.close(fig)


### Fig3a:
def generate_fig3a(ci_accuracy):
	# generate Fig3 from trunc_CI vs est_CI processed matrix
	fig, axs = plt.subplots(1, 1, figsize=(36, 12))
	plt.rcParams.update({'font.size': 40})
	temp_df=ci_accuracy.copy()
	temp_df['name']=['k='+str(x) for x in ci_accuracy.index]
	temp_df.plot(x='name', y=list(ci_accuracy.columns), kind="bar", ax=axs)
	axs.legend(bbox_to_anchor=(1.00, 1.02), title="Depth")
	axs.set(xlabel=None)
	axs.set_ylabel("Consistency ratio") # ,loc="top", rotation=0
	plt.setp(axs.get_xticklabels(), rotation=0)
	fig.savefig("Fig3a_compare_ci.png", dpi=200)
	plt.close(fig)



### Fig3bcd:
def generate_fig3bcd(ref_size_list, trunc_time, est_time):
	fig, axs = plt.subplots(1, 3, figsize=(36, 12))
	plt.rcParams.update({'font.size': 25})
	# space usage, similar to Fig1
	trunc_line = ref_size_list[len(ref_size_list) - 1]
	est_list = list(np.cumsum(ref_size_list))
	input_x = [x + 1 for x in list(range(len(ref_size_list)))]
	sb.lineplot(x=input_x, y=[trunc_line] * len(k_sizes), ax=axs[0], color="orange", linewidth=2, marker="o",
	            markersize=15)
	sb.lineplot(x=input_x, y=est_list, ax=axs[0], color="black", linewidth=2, marker="o", markersize=15)
	legend_element = [
		Line2D([0], [0], marker='o', color='orange', label='CMash', markerfacecolor='orange', markersize=15),
		Line2D([0], [0], marker='o', color='black', label='Classic MinHash', markerfacecolor='black', markersize=15)]
	axs[0].set(ylabel="Database size / MB", xlabel="number of k values")
	# time usage
	input_depth = list(trunc_time['depth'])
	time_trunc_ref = [math.ceil(x / 60) for x in list(trunc_time['ref_time'])]
	time_trunc_run = [math.ceil(x / 60) for x in list(trunc_time['running_time'])]
	time_est_ref = [math.ceil(x / 60) for x in list(est_time['ref_time'])]
	time_est_run = [math.ceil(x / 60) for x in list(est_time['running_time'])]
	### generate plot
	input_x = [x.replace('m', '') for x in input_depth]
	sb.lineplot(x=input_x, y=time_trunc_ref, ax=axs[1], color="orange", linewidth=2, marker="o", markersize=15)
	sb.lineplot(x=input_x, y=time_est_ref, ax=axs[1], color="black", linewidth=2, marker="o", markersize=15)
	axs[1].legend(handles=legend_element, bbox_to_anchor=(0.75, 1.15), fontsize=22)
	axs[1].set(ylabel="CPU time / min", xlabel="Metagenomic depth / Million reads")
	axs[1].set_ylim([0, 2000])
	sb.lineplot(x=input_x, y=time_trunc_run, ax=axs[2], color="orange", linewidth=2, marker="o", markersize=15)
	sb.lineplot(x=input_x, y=time_est_run, ax=axs[2], color="black", linewidth=2, marker="o", markersize=15)
	axs[2].set(ylabel="CPU time / min", xlabel="Metagenomic depth / Million reads")
	axs[2].set_ylim([0, 2000])
	
	fig.savefig("Fig3_bcd_compare_time_space.png", dpi=200)
	plt.close(fig)



# local test
def local_variable():
	# activate local variables for test purpose
	target_dir = "/Users/shaopeng/Desktop/final_output"
	os.chdir(target_dir)
	os.listdir()



# generate plots
if __name__ == '__main__':
	print("Generating figures in the manuscript")
	
	# workdir
	out_dir = os.getcwd()
	f2_dir = out_dir+"/fig_1f_2_input"
	f3_dir = out_dir+"/fig_1f_3_input"
	
	# for fig2: read all csv files
	os.chdir(f2_dir)
	k_sizes = [15, 20, 25, 30, 35, 40, 45, 50, 55, 60]  # for fig2
	dict_est = read_all_k_df('est', k_sizes)
	dict_gt = read_all_k_df('gt', k_sizes)
	dict_trunc = read_all_k_df('trunc', k_sizes)
	### for fig2a
	df_gt = dict_to_df_by_k(dict_gt, k_sizes)
	### for fig2b
	df_trunc = dict_to_df_by_k(dict_trunc, k_sizes)
	
	# for fig3: read all csv/txt files
	os.chdir(f3_dir)
	filename_key=[str(x) + 'm' for x in list(range(2, 21, 2))] #those BBMap files
	meta_ksizes = [20, 25, 30, 35, 40, 45, 50, 55, 60]  # for fig3
	ci_accuracy = collect_ci_dif_data(k_range=meta_ksizes, name_key=filename_key, cutoff=0.05)
	### the cumsum of ref size
	ref_size_list = list(pd.read_table("cum_space_est_CI.txt", header=None).iloc[:,0]) #total trunc size is the last element in est_CI ref size
	ref_size_list = [math.ceil(x/1024/1024) for x in ref_size_list] #change to MB
	### time record files
	trunc_time = pd.read_table("time_summary_trunc_CI.txt", header=0)
	est_time = pd.read_table("time_summary_est_CI.txt", header=0)
	
	# ready for plot
	os.chdir(out_dir)
	### fig1f.left
	[abs_dif, rela_dif] = get_dif_matrix_and_fig1f(dict_trunc, dict_gt, k_sizes)
	### fig1f.right
	generate_fig1f_right(ref_size_list=ref_size_list, trunc_time=trunc_time, est_time=est_time)
	### fig2
	generate_fig2(df_trunc=df_trunc, df_gt=df_gt, df_abs_dif=abs_dif, df_rela_dif=rela_dif)
	### fig3a
	generate_fig3a(ci_accuracy=ci_accuracy)
	### fig3bcd
	generate_fig3bcd(ref_size_list=ref_size_list, trunc_time=trunc_time, est_time=est_time)
	
	
	print("Ploting finished")
	
	
	
	
	
	
	


