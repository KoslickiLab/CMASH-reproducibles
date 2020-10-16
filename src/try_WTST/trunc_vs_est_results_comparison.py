# need to connect to previous code for a whole process

import argparse
import os
import glob
import pandas as pd
import numpy as np
import sys
import seaborn as sb
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

# temp variables
input_range="10-60-5"

parser = argparse.ArgumentParser(description='Generating plots')
parser.add_argument('-n', '--name', type=str, help="name for plot label", default="input data")
args = parser.parse_args()
name = args.name


start=input_range.split("-")[0]
end=input_range.split("-")[1]
gap=input_range.split("-")[2]

#note: df2 = df2.where(np.triu(np.ones(df2.shape), 1).astype(np.bool))
def upper_triangle(input_df):
	out_df = input_df.where(np.triu(np.ones(input_df.shape), 1).astype(np.bool))
	return out_df

# generate a list of variance number / percentage for 3 matched files
def generate_list(trunc, est, gt, symmatric=False):
	df_trunc = pd.read_csv(trunc, header=0, index_col=0)
	df_est = pd.read_csv(est, header=0, index_col=0)
	df_gt = pd.read_csv(gt, header=0, index_col=0)  #used as base
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

for k in range(int(start), int(end), int(gap)):
	temp_colname = 'k='+str(k)
	v_trunc, v_est, v_gt, vd_trunc_gt, vd_est_gt, vd_trunc_est, vr_trunc_gt, vr_est_gt, vr_trunc_est = generate_list("trunc_WJI_k"+str(k)+".csv", "est_WJI_k"+str(k)+".csv", "GT_WJI_k"+str(k)+".csv" ,symmatric=True)
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
print("Generating plot for "+name)

# boxplot
fig1, axs = plt.subplots(3)
sb.boxplot(data=p_trunc, color="red", ax=axs[0])
sb.boxplot(data=p_est, color="cyan", ax=axs[1])
sb.boxplot(data=p_gt, color="yellow", ax=axs[2])
axs[0].set(ylabel="Trunc_WJI")
axs[1].set(ylabel="Est_WJI")
axs[2].set(ylabel="GroundTruth_WJI")
fig1.tight_layout(rect=[0, 0.03, 1, 0.9])
fig1.suptitle("Boxplot of 3 WJI")
fig1.savefig("Boxplot_of_" + name + ".png", dpi=200)
plt.close(fig1)

### merged swam plot
fig3, axs3 = plt.subplots()
sb.swarmplot(data=p_trunc, color="red", size=3, linewidth=0.3, ax=axs3)
sb.swarmplot(data=p_est, color="cyan", size=3, linewidth=0.3, ax=axs3)
sb.swarmplot(data=p_gt, color="black", size=3, linewidth=0.3, ax=axs3)
axs3.set(ylabel="WJI value")
legend_element = [Line2D([0], [0], marker='o', color='g', label='Trunc', markerfacecolor='red', markersize=5),
				  Line2D([0], [0], marker='o', color='g', label='Est', markerfacecolor='cyan', markersize=5),
				  Line2D([0], [0], marker='o', color='g', label='GT', markerfacecolor='black', markersize=5)]
axs3.legend(handles=legend_element, loc='best')
fig3.suptitle("Merged Swarmplot")
fig3.savefig("Merged_swarmplot_of_" + name + ".png", dpi=200)
plt.close(fig3)



# dif plot, temp_gt.min = 0.00003
fig4, axs4 = plt.subplots(3, 2, figsize=(16,12))
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
	medians=df.median()
	vertical_offset = max(medians.median(), 0.005)
	for xtick in ax.get_xticks():
		ax.text(xtick, medians[xtick] + vertical_offset, round(medians[xtick], 3), horizontalalignment='center',
				size='x-small', color='navy', weight='semibold')

# boxplot of WJI change: truncated vs groundtruth
fig5, axs = plt.subplots(3, figsize=(9, 9))
fig5.suptitle("Boxplot of WJI dif")
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
fig5.savefig("Boxplot_of_CI_dif_" + name + ".png", dpi=300)
plt.close(fig5)