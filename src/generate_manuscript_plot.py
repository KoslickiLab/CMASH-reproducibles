import os
import pandas as pd
import numpy as np
import seaborn as sb
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import math
import statistics
import string



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
def updated_ci_df_for_box_3a(k_range, name_key):
	# this will return a dict, key is the depth 2m/4m etc, each value is a df that stores the absolute difference
	out_dict = dict()
	for key in name_key:
		df_key = pd.DataFrame()
		f_name = "trunc_CI_results_" + key + ".csv"
		# truncated CI file contains estimation for all k values (simoutaneously)
		trunc_df = pd.read_csv(f_name, header=0, index_col=0)
		for k in k_range:
			# now is the est_CI by classic minhash
			f_name = "est_CI_results_" + key + "_k" + str(k) + ".csv"
			temp_k_df = pd.read_csv(f_name, header=0, index_col=0)
			temp_trunc_df = trunc_df[["k=" + str(k)]]
			# remember to merge by colname: 2 files have different sort!!!
			temp_check = pd.concat([temp_k_df, temp_trunc_df], axis=1, join="inner")
			df_key['k='+str(k)]=abs(temp_check.iloc[:, 0] - temp_check.iloc[:, 1])
		# put depth_df to out_dict
		out_dict[key] = df_key
	
	return out_dict
	

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
	plt.rcParams.update({'font.size': 12})
	# ab: swarmplot of GT and Trunc JI
	sb.swarmplot(data=df_gt, ax=axs[0,0], color="black", size=1, linewidth=0.3)
	axs[0,0].set(ylabel="Ground truth Jaccard indices")
	axs[0,0].set_ylim([0, 1])
	sb.swarmplot(data=df_trunc, ax=axs[0, 1], color="black", size=1, linewidth=0.3)
	axs[0, 1].set(ylabel="CMash Jaccard indices")
	axs[0, 1].set_ylim([0, 1])
	# cd: absolute error and relative error
	sb.boxplot(data=df_abs_dif, ax=axs[1,0], color="black")
	axs[1,0].set(ylabel="CMash JI minus Ground truth JI")
	axs[1,0].set_ylim([-0.3, 0.3])
	sb.boxplot(data=df_rela_dif, ax=axs[1,1], color="black")
	axs[1,1].set(ylabel="Relative error: CMash JI to Ground truth JI")
	axs[1,1].set_ylim([-0.3, 0.3])
	for i in [0,1]:
		for j in [0,1]:
			plt.setp(axs[i,j].get_xticklabels(), rotation=45)
	#fig.suptitle(plot_name)
	fig.savefig("Fig2_compare_trunc_gt_ji.png", dpi=200)
	plt.close(fig)

def horizontal_fig2(df_trunc, df_gt, df_abs_dif, df_rela_dif, font=30, dpi=100):
	# generate Fig2 from Trunc/GT processed matrix
	fig, axs = plt.subplots(1,4, figsize=(40, 10))
	plt.rcParams.update({'font.size': font})
	# ab: swarmplot of GT and Trunc JI
	sb.swarmplot(data=df_gt, ax=axs[0], color="black", size=1, linewidth=0.3)
	axs[0].set(ylabel="Ground truth Jaccard indices")
	axs[0].set_ylim([0, 1])
	sb.swarmplot(data=df_trunc, ax=axs[1], color="black", size=1, linewidth=0.3)
	axs[1].set(ylabel="CMash Jaccard indices")
	axs[1].set_ylim([0, 1])
	# cd: absolute error and relative error
	sb.boxplot(data=df_abs_dif, ax=axs[2], color="black")
	axs[2].set(ylabel="CMash JI minus Ground truth JI")
	axs[2].set_ylim([-0.3, 0.3])
	sb.boxplot(data=df_rela_dif, ax=axs[3], color="black")
	axs[3].set(ylabel="Relative error: CMash JI to Ground truth JI")
	axs[3].set_ylim([-0.3, 0.3])
	letter_mark=['(a)','(b)','(c)','(d)']
	for i in [0,1,2,3]:
		plt.setp(axs[i].get_xticklabels(), rotation=45)
		# add letter mark
		axs[i].text(-0.1, 1.1, letter_mark[i], transform=axs[i].transAxes, size=20, weight='bold')
	#fig.suptitle(plot_name)
	fig.tight_layout()
	fig.savefig("f2_horizontal.png", dpi=dpi)
	plt.close(fig)


def ismb_fig2(df_trunc, df_est, df_gt, df_trunc_abs_dif, df_trunc_rela_dif, df_est_abs_dif, df_est_rela_dif, font=30,
              dpi=100):
	# will put 3 separate figures first
	# generate Fig2 from Est(MH)/Trunc(CMash)/GT matrix
	# 2 figures: 1 + 4
	plt.rcParams.update({'font.size': font})
	
	### part1: swarmplot GT only
	fig, axs = plt.subplots(figsize=(10, 10))
	sb.swarmplot(data=df_gt, ax=axs, color="black", size=1, linewidth=0.3)
	axs.set(ylabel="Ground truth Jaccard indices")
	axs.set_ylim([0, 1])
	plt.setp(axs.get_xticklabels(), rotation=45)
	# add letter mark
	axs.text(-0.1, 1.1, '(a)', transform=axs.transAxes, size=20, weight='bold')
	fig.tight_layout()
	fig.savefig("f2_part1.png", dpi=dpi / 2)
	plt.close(fig)
	
	### part2: 2x2, CMash vs GT, Est vs GT
	fig, axs = plt.subplots(2, 2, figsize=(20, 20))
	sb.boxplot(data=df_trunc_abs_dif, ax=axs[0, 0], color="black")
	axs[0, 0].set(ylabel="CMash JI minus Ground truth JI")
	axs[0, 0].set_ylim([-0.3, 0.3])
	sb.boxplot(data=df_trunc_rela_dif, ax=axs[0, 1], color="black")
	axs[0, 1].set(ylabel="Relative error: CMash JI to Ground truth JI")
	axs[0, 1].set_ylim([-0.3, 0.3])
	sb.boxplot(data=df_est_abs_dif, ax=axs[1, 0], color="black")
	axs[1, 0].set(ylabel="MinHash JI minus Ground truth JI")
	axs[1, 0].set_ylim([-0.3, 0.3])
	sb.boxplot(data=df_est_rela_dif, ax=axs[1, 1], color="black")
	axs[1, 1].set(ylabel="Relative error: MinHash JI to Ground truth JI")
	axs[1, 1].set_ylim([-0.3, 0.3])
	### letter mark
	letter_mark = [['(b)', '(c)'], ['(d)', '(e)']]
	for i in [0, 1]:
		for j in [0, 1]:
			plt.setp(axs[i, j].get_xticklabels(), rotation=45)
			axs[i, j].text(-0.1, 1.1, letter_mark[i][j], transform=axs[i, j].transAxes, size=20, weight='bold')
	fig.tight_layout()
	fig.savefig("f2_part2.png", dpi=dpi)
	plt.close(fig)


### Fig3a
def updated_fig3a_box(f3a_dict, dict_key, drop_last=True, font=25, dpi=100):
	# generate box plot for a selected df (depth) in f3a_dict
	fig, axs = plt.subplots(1, 1, figsize=(36, 12))
	plt.rcParams.update({'font.size': font})
	temp_df = f3a_dict[dict_key]
	# last column is untruncated, do dropped
	if drop_last:
		temp_df=temp_df.iloc[: , :-1]
	# generate box plot or swarmplot
	sb.swarmplot(data=temp_df ,ax=axs, size=3, color="grey")
	sb.boxplot(data=temp_df, ax=axs, color='bisque', showfliers = False)
	axs.set_ylabel("Absolute difference in CI values")  # ,loc="top", rotation=0
	axs.set(title='Boxplot for CMash CI vs MinHash CI')
	axs.set_ylim([0, 0.08])
	plt.setp(axs.get_xticklabels(), rotation=0)
	fig.savefig("Fig3a_compare_ci.png", dpi=dpi)
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
	target_dir = "/Users/shaopeng/OneDrive - The Pennsylvania State University/PSU_academic/Koslicki_lab/publication/conference/WABI_2021"
	os.chdir(target_dir)
	os.listdir()
	# refresh fig2
	horizontal_fig2(df_trunc=df_trunc, df_gt=df_gt, df_abs_dif=abs_dif, df_rela_dif=rela_dif, font=30, dpi=100)
	# refresh fig3a
	updated_fig3a_box(f3a_dict=f3a_dict, dict_key='10m', font=30, dpi=200)
	
def code_storage():
	### for CMash vs GT, Est vs GT box plot:
	fig, axs = plt.subplots(1, 1, figsize=(36, 12))
	temp_df = df_trunc_abs_dif
	sb.swarmplot(data=temp_df, ax=axs, size=3, color="grey")
	sb.boxplot(data=temp_df, ax=axs, color='bisque', showfliers=False)
	axs.set_ylabel("Difference in JI values")  # ,loc="top", rotation=0
	axs.set(title='Boxplot for CMash JI vs GroundTruth JI')
	axs.set_ylim([-0.05, 0.05])
	plt.setp(axs.get_xticklabels(), rotation=0)
	fig.savefig("Fig2_add_ci.png", dpi=dpi)
	plt.close(fig)
	
	fig, axs = plt.subplots(1, 1, figsize=(36, 12))
	temp_df = df_est_abs_dif
	sb.swarmplot(data=temp_df, ax=axs, size=3, color="grey")
	sb.boxplot(data=temp_df, ax=axs, color='bisque', showfliers=False)
	axs.set_ylabel("Difference in JI values")  # ,loc="top", rotation=0
	axs.set(title='Boxplot for MinHash JI vs GroundTruth JI')
	axs.set_ylim([-0.05, 0.05])
	plt.setp(axs.get_xticklabels(), rotation=0)
	fig.savefig("Fig2_add_mh.png", dpi=dpi)
	plt.close(fig)


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
	### Classic MH
	df_est = dict_to_df_by_k(dict_est, k_sizes)
	
	# for fig3: read all csv/txt files
	os.chdir(f3_dir)
	filename_key=[str(x) + 'm' for x in list(range(2, 21, 2))] #those BBMap files
	meta_ksizes = [20, 25, 30, 35, 40, 45, 50, 55, 60]  # for fig3
	
	# fig3a box data
	f3a_dict = updated_ci_df_for_box_3a(k_range=meta_ksizes, name_key=filename_key)
	
	
	### the cumsum of ref size
	ref_size_list = list(pd.read_table("cum_space_est_CI.txt", header=None).iloc[:,0]) #total trunc size is the last element in est_CI ref size
	ref_size_list = [math.ceil(x/1024/1024) for x in ref_size_list] #change to MB
	### time record files
	trunc_time = pd.read_table("time_summary_trunc_CI.txt", header=0)
	est_time = pd.read_table("time_summary_est_CI.txt", header=0)
	
	# ready for plot
	os.chdir(out_dir)
	### fig1f.left
	### for CMash vs GT
	[abs_dif, rela_dif] = get_dif_matrix_and_fig1f(dict_trunc, dict_gt, k_sizes)
	
	### for Est vs GT
	[est_abs_dif, est_rela_dif] = get_dif_matrix_and_fig1f(dict_est, dict_gt, k_sizes)
	
	
	### fig1f.right
	generate_fig1f_right(ref_size_list=ref_size_list, trunc_time=trunc_time, est_time=est_time)
	### fig2
	#generate_fig2(df_trunc=df_trunc, df_gt=df_gt, df_abs_dif=abs_dif, df_rela_dif=rela_dif)
	#horizontal_fig2(df_trunc=df_trunc, df_gt=df_gt, df_abs_dif=abs_dif, df_rela_dif=rela_dif)
	ismb_fig2(df_trunc=df_trunc, df_est=df_est, df_gt=df_gt, df_trunc_abs_dif=abs_dif, df_trunc_rela_dif=rela_dif,
	          df_est_abs_dif=est_abs_dif, df_est_rela_dif=est_rela_dif)
	### fig3a
	updated_fig3a_box(f3a_dict=f3a_dict, dict_key='10m')
	
	### fig3bcd
	generate_fig3bcd(ref_size_list=ref_size_list, trunc_time=trunc_time, est_time=est_time)
	
	
	print("Ploting finished")
	
	



	




	
	
	


	
	