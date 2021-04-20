import argparse
import os
import glob
import pandas as pd
import numpy as np
import sys
import seaborn as sb
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import re

parser = argparse.ArgumentParser(description='Generating plots')
parser.add_argument('-n', '--name', type=str, help="name for plot label", default="input_data")
parser.add_argument('-t', '--tril', type=str, help="Use lower triangle of an matrix", default="False")
parser.add_argument('-g', '--k_range', type=str, help="k-mer range", default="10-60-5")
parser.add_argument('-k', '--maxk', type=int, help="k-mer range", default=60)
args = parser.parse_args()
name = args.name
input_range = args.k_range
maxk = args.maxk
use_tril = args.tril


# pipe start
def parsenumlist(k_sizes_str: str):
	"""
	Parses a string like 10-21-1 and turn it into a list like [10, 11, 12,...,21]
	:param k_sizes_str: the <start>-<end>-<increment> string
	:type k_sizes_str: str
	:return: list of k-mer sizes
	:rtype: list
	"""
	m = re.match(r'(\d+)(?:-(\d+))?(?:-(\d+))?$', k_sizes_str)
	# ^ (or use .split('-'). anyway you like.)
	if not m:
		raise argparse.ArgumentTypeError(
			"'" + k_sizes_str + "' is not a range of number. Expected forms like '1-5' or '2' or '10-15-2'.")
	start = int(m.group(1))
	end = int(m.group(2))
	if m.group(3):
		increment = int(m.group(3))
	else:
		increment = 1
	return list(range(start, end + 1, increment))


k_sizes = parsenumlist(input_range)
use_tril = use_tril == "True"


# clean matrix for symmatric input
def upper_triangle(input_df):
	if use_tril:
		print("Ploting by the lower triangle of the input matrix!!!")
		out_df = input_df.where(np.tril(np.ones(input_df.shape), -1).astype(np.bool))
	else:
		out_df = input_df.where(np.triu(np.ones(input_df.shape), 1).astype(np.bool))
	return out_df


# plot est / trunc / gt comparison
# generate a list of variance number / percentage for 3 matched files
def generate_list(trunc, est, gt, symmatric=False):
	df_trunc = pd.read_csv(trunc, header=0, index_col=0)
	df_est = pd.read_csv(est, header=0, index_col=0)
	df_gt = pd.read_csv(gt, header=0, index_col=0)  # used as base
	dif_trunc_gt = df_trunc - df_gt
	dif_est_gt = df_est - df_gt
	dif_trunc_est = df_trunc - df_est
	ratio_trunc_gt = dif_trunc_gt / (df_gt + 0.0001)
	ratio_est_gt = dif_est_gt / (df_gt + 0.0001)
	ratio_trunc_est = dif_trunc_est / (df_est + 0.0001)
	
	# if the matrix is self-matched (symmatric, diagnol=1), remove those unuseful records
	if symmatric:
		df_trunc = upper_triangle(df_trunc)
		df_est = upper_triangle(df_est)
		df_gt = upper_triangle(df_gt)
		dif_trunc_gt = upper_triangle(dif_trunc_gt)
		dif_est_gt = upper_triangle(dif_est_gt)
		dif_trunc_est = upper_triangle(dif_trunc_est)
		ratio_trunc_gt = upper_triangle(ratio_trunc_gt)
		ratio_est_gt = upper_triangle(ratio_est_gt)
		ratio_trunc_est = upper_triangle(ratio_trunc_est)
	
	value_trunc = df_trunc.stack().to_numpy()
	value_est = df_est.stack().to_numpy()
	value_gt = df_gt.stack().to_numpy()
	value_dif_trunc_gt = dif_trunc_gt.stack().to_numpy()
	value_dif_est_gt = dif_est_gt.stack().to_numpy()
	value_dif_trunc_est = dif_trunc_est.stack().to_numpy()
	value_ratio_trunc_gt = ratio_trunc_gt.stack().to_numpy()
	value_ratio_est_gt = ratio_est_gt.stack().to_numpy()
	value_ratio_trunc_est = ratio_trunc_est.stack().to_numpy()
	
	return value_trunc, value_est, value_gt, value_dif_trunc_gt, value_dif_est_gt, value_dif_trunc_est, value_ratio_trunc_gt, value_ratio_est_gt, value_ratio_trunc_est


# read csv and merge the results into dataframe for plotting
# p: ready to plot
p_trunc = pd.DataFrame()
p_est = pd.DataFrame()
p_gt = pd.DataFrame()
pd_trunc_gt = pd.DataFrame()
pd_est_gt = pd.DataFrame()
pd_trunc_est = pd.DataFrame()
pr_trunc_gt = pd.DataFrame()
pr_est_gt = pd.DataFrame()
pr_trunc_est = pd.DataFrame()

for k in k_sizes:
	temp_colname = 'k=' + str(k)
	v_trunc, v_est, v_gt, vd_trunc_gt, vd_est_gt, vd_trunc_est, vr_trunc_gt, vr_est_gt, vr_trunc_est = generate_list(
		"trunc_JI_k" + str(k) + ".csv", "est_JI_k" + str(k) + ".csv", "GroundTruth_JI_k" + str(k) + ".csv",
		symmatric=True)
	p_trunc[temp_colname] = v_trunc
	p_est[temp_colname] = v_est
	p_gt[temp_colname] = v_gt
	pd_trunc_gt[temp_colname] = vd_trunc_gt
	pd_est_gt[temp_colname] = vd_est_gt
	pd_trunc_est[temp_colname] = vd_trunc_est
	pr_trunc_gt[temp_colname] = vr_trunc_gt
	pr_est_gt[temp_colname] = vr_est_gt
	pr_trunc_est[temp_colname] = vr_trunc_est

# generate plots
print("Generating plot for " + name)

# boxplot
# fig1, axs = plt.subplots(3)
# sb.boxplot(data=p_trunc, color="red", ax=axs[0])
# sb.boxplot(data=p_est, color="cyan", ax=axs[1])
# sb.boxplot(data=p_gt, color="yellow", ax=axs[2])
# axs[0].set(ylabel="Trunc_JI")
# axs[1].set(ylabel="Est_JI")
# axs[2].set(ylabel="GroundTruth_JI")
# fig1.tight_layout(rect=[0, 0.03, 1, 0.9])
# fig1.suptitle("Boxplot of 3 JI")
# fig1.savefig("Boxplot_of_" + name + ".png", dpi=200)
# plt.close(fig1)

### merged swam plot
fig3, axs3 = plt.subplots()
sb.swarmplot(data=p_trunc, color="red", size=3, linewidth=0.3, ax=axs3)
sb.swarmplot(data=p_est, color="cyan", size=3, linewidth=0.3, ax=axs3)
sb.swarmplot(data=p_gt, color="black", size=3, linewidth=0.3, ax=axs3)
axs3.set(ylabel="JI value")
legend_element = [Line2D([0], [0], marker='o', color='g', label='Trunc', markerfacecolor='red', markersize=5),
                  Line2D([0], [0], marker='o', color='g', label='Est', markerfacecolor='cyan', markersize=5),
                  Line2D([0], [0], marker='o', color='g', label='GT', markerfacecolor='black', markersize=5)]
axs3.legend(handles=legend_element, loc='best')
fig3.suptitle("Merged Swarmplot")
fig3.savefig("Merged_swarmplot_of_" + name + ".png", dpi=200)
plt.close(fig3)

# dif plot, temp_gt.min = 0.00003
fig4, axs4 = plt.subplots(3, 2, figsize=(16, 12))
sb.swarmplot(data=pd_trunc_gt, ax=axs4[0, 0], size=3, linewidth=0.3)
sb.swarmplot(data=pr_trunc_gt, ax=axs4[0, 1], size=3, linewidth=0.3)
sb.swarmplot(data=pd_est_gt, ax=axs4[1, 0], size=3, linewidth=0.3)
sb.swarmplot(data=pr_est_gt, ax=axs4[1, 1], size=3, linewidth=0.3)
sb.swarmplot(data=pd_trunc_est, ax=axs4[2, 0], size=3, linewidth=0.3)
sb.swarmplot(data=pr_trunc_est, ax=axs4[2, 1], size=3, linewidth=0.3)
axs4[0, 0].set(ylabel="Truncated vs GroundTruth")
axs4[1, 0].set(ylabel="Estimated vs GroundTruth")
axs4[2, 0].set(ylabel="Truncated vs Estimated")
axs4[0, 0].title.set_text("Change by number")
axs4[0, 1].title.set_text("Change by ratio")
fig4.suptitle("Pairwise comparison")
fig4.savefig("Pairwise_comparison_of_" + name + ".png", dpi=200)
plt.close(fig4)


def label_median(df, ax):
	medians = df.median()
	vertical_offset = max(medians.median(), 0.005)
	for xtick in ax.get_xticks():
		ax.text(xtick, medians[xtick] + vertical_offset, round(medians[xtick], 3), horizontalalignment='center',
		        size='x-small', color='navy', weight='semibold')


# boxplot of JI change: truncated vs groundtruth
fig5, axs = plt.subplots(3, figsize=(9, 9))
fig5.suptitle("Boxplot of JI dif")
sb.boxplot(data=pd_trunc_gt, color="red", ax=axs[0])
label_median(pd_trunc_gt, axs[0])
sb.boxplot(data=pd_est_gt, color="cyan", ax=axs[1])
label_median(pd_est_gt, axs[1])
sb.boxplot(data=pd_trunc_est, color="orange", ax=axs[2])
label_median(pd_trunc_est, axs[2])
axs[0].set(ylabel="Truncated-Groundtruth")
axs[1].set(ylabel="Estimated-Groundtruth")
axs[2].set(ylabel="Truncated-Estimated")
fig5.tight_layout(rect=[0, 0.03, 1, 0.9])
fig5.savefig("Boxplot_of_JI_dif_" + name + ".png", dpi=300)
plt.close(fig5)

# uncorrected dif plot
# dif plot, temp_gt.min = 0.00003
fig4, axs4 = plt.subplots(2, 2, figsize=(16, 12))
sb.swarmplot(data=pd_trunc_gt, ax=axs4[0, 0], size=3, linewidth=0.3)
sb.swarmplot(data=pr_trunc_gt, ax=axs4[0, 1], size=3, linewidth=0.3)
sb.boxplot(data=pd_trunc_gt, ax=axs4[1, 0], color="red")
label_median(pd_trunc_gt, axs4[1, 0])
sb.boxplot(data=pr_trunc_gt, ax=axs4[1, 1], color="cyan")
label_median(pr_trunc_gt, axs4[1, 1])
axs4[0, 0].set(ylabel="Trunc - GT")
axs4[0, 0].title.set_text("Change by number")
axs4[0, 1].title.set_text("Change by ratio")
axs4[1, 0].set(ylabel="Trunc - GT")
axs4[1, 0].title.set_text("Change by number")
axs4[1, 1].title.set_text("Change by ratio")
fig4.suptitle("Deviance_before_bias_correcton")
fig4.savefig("Plot_of_uncorrected_trunc_" + name + ".png", dpi=200)
plt.close(fig4)


def generate_bias_list(trunc, bias, gt, symmatric=False):
	df_trunc = pd.read_csv(trunc, header=0, index_col=0)
	df_bias = pd.read_csv(bias, header=0, index_col=0)
	df_gt = pd.read_csv(gt, header=0, index_col=0)  # used as base
	df_multiple = df_gt * df_bias
	df_dif = df_trunc - df_multiple
	ratio_dif_gt = df_dif / (df_gt + 0.0001)
	
	# if the matrix is self-matched (symmatric, diagnol=1), remove those unuseful records
	if symmatric:
		df_trunc = upper_triangle(df_trunc)
		df_gt = upper_triangle(df_gt)
		df_dif = upper_triangle(df_dif)
		df_multiple = upper_triangle(df_multiple)
		ratio_dif_gt = upper_triangle(ratio_dif_gt)
	
	value_trunc = df_trunc.stack().to_numpy()
	value_gt = df_gt.stack().to_numpy()
	value_dif = df_dif.stack().to_numpy()
	value_ratio_dif_gt = ratio_dif_gt.stack().to_numpy()
	df_multiple = df_multiple.stack().to_numpy()
	
	return value_trunc, value_gt, value_dif, value_ratio_dif_gt, df_multiple


p_trunc = pd.DataFrame()
p_gt = pd.DataFrame()
p_dif = pd.DataFrame()
pr_dif = pd.DataFrame()
p_multiple = pd.DataFrame()

if maxk in k_sizes:
	k_sizes.remove(maxk)

for k in k_sizes:
	temp_colname = 'k=' + str(k)
	v_trunc, v_gt, v_dif, vr_dif, v_multiple = generate_bias_list("trunc_JI_k" + str(k) + ".csv",
	                                                              "bias_factor_k" + str(k) + "_to_k" + str(
		                                                              maxk) + ".csv",
	                                                              "GroundTruth_JI_k" + str(k) + ".csv", symmatric=True)
	p_trunc[temp_colname] = v_trunc
	p_gt[temp_colname] = v_gt
	p_dif[temp_colname] = v_dif
	pr_dif[temp_colname] = vr_dif
	p_multiple[temp_colname] = v_multiple

# plot the bias factor
### merged swam plot
fig3, axs3 = plt.subplots()
sb.swarmplot(data=p_trunc, color="red", size=3, linewidth=0.3, ax=axs3)
sb.swarmplot(data=p_multiple, color="cyan", size=3, linewidth=0.3, ax=axs3)
sb.swarmplot(data=p_gt, color="black", size=3, linewidth=0.3, ax=axs3)
axs3.set(ylabel="JI value")
legend_element = [Line2D([0], [0], marker='o', color='g', label='Trunc', markerfacecolor='red', markersize=5),
                  Line2D([0], [0], marker='o', color='g', label='GT_times_bias', markerfacecolor='cyan', markersize=5),
                  Line2D([0], [0], marker='o', color='g', label='GT', markerfacecolor='black', markersize=5)]
axs3.legend(handles=legend_element, loc='best')
fig3.suptitle("Merged Swarmplot")
fig3.savefig("Merged_swarmplot_of_GT_trunc_bias_" + name + ".png", dpi=200)
plt.close(fig3)

# dif plot, temp_gt.min = 0.00003
fig4, axs4 = plt.subplots(2, 2, figsize=(16, 12))
sb.swarmplot(data=p_dif, ax=axs4[0, 0], size=3, linewidth=0.3)
sb.swarmplot(data=pr_dif, ax=axs4[0, 1], size=3, linewidth=0.3)
sb.boxplot(data=p_dif, ax=axs4[1, 0], color="red")
label_median(p_dif, axs4[1, 0])
sb.boxplot(data=pr_dif, ax=axs4[1, 1], color="cyan")
label_median(p_dif, axs4[1, 1])
axs4[0, 0].set(ylabel="Trunc - GT*bias")
axs4[0, 0].title.set_text("Change by number")
axs4[0, 1].title.set_text("Change by ratio")
axs4[1, 0].set(ylabel="Trunc - GT*bias")
axs4[1, 0].title.set_text("Change by number")
axs4[1, 1].title.set_text("Change by ratio")
fig4.suptitle("Deviance_after_bias_correcton")
fig4.savefig("Plot_of_corrected_trunc_" + name + ".png", dpi=200)
plt.close(fig4)

### plot the asymmatric matrix:
dif_matrix = pd.DataFrame()
rdif_matrix = pd.DataFrame()
bias_matrix = pd.DataFrame()

for k in k_sizes:
	temp_colname = 'k=' + str(k)
	df = pd.read_csv(f"trunc_JI_k{k}.csv", header=0, index_col=0)
	dim = df.shape[0]
	# get all pairwise difference: upper - lower
	dif_list = []
	rdif_list = []
	for i in range(dim):
		for j in range(dim):
			if i < j:
				temp = df.iloc[i, j] - df.iloc[j, i]
				rtemp = temp / (df.iloc[j, i] + 0.00001)
				dif_list.append(temp)
				rdif_list.append(rtemp)
	dif_matrix[temp_colname] = dif_list
	rdif_matrix[temp_colname] = rdif_list
	# add bias matrix
	bias = pd.read_csv("bias_factor_k" + str(k) + "_to_k" + str(maxk) + ".csv", header=0, index_col=0)
	bias = upper_triangle(bias)
	bias = bias - 1
	bias_list = bias.stack().to_numpy()
	bias_matrix[temp_colname] = bias_list

### plot
fig4, axs4 = plt.subplots(2, 2, figsize=(16, 12))
sb.boxplot(data=dif_matrix, ax=axs4[0, 0], color="red")
label_median(dif_matrix, axs4[0, 0])
sb.boxplot(data=bias_matrix, ax=axs4[1, 1], color="cyan")
label_median(bias_matrix, axs4[1, 1])
sb.boxplot(data=rdif_matrix, ax=axs4[1, 0], color="cyan")
label_median(rdif_matrix, axs4[1, 0])
axs4[0, 0].set(ylabel="Dif by value")
axs4[0, 0].title.set_text("Trunc_JI matrix dif by value")
axs4[1, 1].title.set_text("Bias factor minus 1")
axs4[1, 0].title.set_text("Trunc_JI matrix dif by ratio")
axs4[1, 0].set(ylabel="Dif by ratio")
fig4.suptitle("Matrix difference and bias factor")
fig4.savefig("Plot_of_matrix_dif_of_trunc_JI_" + name + ".png", dpi=200)
plt.close(fig4)

