#! /usr/bin/env python
import argparse
import os
import glob
import pandas as pd
import numpy as np
import sys
import seaborn as sb
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

### import parameters
parser = argparse.ArgumentParser(description="Summarize the results of Groundtruth_CI, estimated_CI and truncated_CI.")
parser.add_argument('range', type=str, help="Range of k-mer sizes in the formate <start>-<end>-<increment>.")
parser.add_argument('query', type=str, help="Query file of analysis shown as absolute path.")

args = parser.parse_args()
input_range= args.range
start=input_range.split("-")[0]
end=input_range.split("-")[1]
gap=input_range.split("-")[2]
query_file= args.query

### merge all results by each metagenome input
f=open(query_file, 'r')
query_list=[ x.strip() for x in list(f)]
f.close()

### summarize groundtruth_CI
os.chdir("ground_truth_CI")
for line in query_list:
    name=line.split("/")[-1]
    print("Processing ground truth of file: " + name)
    merged_gt_CI=pd.concat([ pd.read_csv(x, header=0, index_col=0) for x in glob.glob("ground_truth_CI_k*"+name+"_results.csv") ] , axis=1, join="inner", sort=True)
    merged_gt_CI=merged_gt_CI[sorted(merged_gt_CI.columns, key=lambda x: int(x[2:]))]
    merged_gt_CI.to_csv("merged_Ground_truth_CI_"+name+".csv", index=True)
    
os.popen("mv merged_Ground_truth_CI_* ../summary")
os.chdir("../estimated_CI/")

### summary estimated_CI
for line in query_list:
    name=line.split("/")[-1]
    print("Processing estimated CI of file: " + name)
    merged_est_CI=pd.concat([ pd.read_csv(x, header=0, index_col=0) for x in glob.glob("estimated_CI_k*"+name+"_results.csv")], axis=1, join="inner", sort=True)
    merged_est_CI=merged_est_CI[sorted(merged_est_CI.columns, key=lambda x: int(x[2:]))]
    merged_est_CI.to_csv("merged_Estimated_CI_"+name+".csv", index=True)

os.popen("mv merged_Estimated_CI_* ../summary")

### the truncated CI was generated in previous running, but need to sort value to make sure they are consistent with the 2 above
os.chdir("../truncated_CI")
for line in query_list:
    name=line.split("/")[-1]
    print("Sorting truncated CI of file: " + name)
    merged_truncat_CI=pd.read_csv("truncation_"+name+"_results.csv", header=0, index_col=0)
    merged_truncat_CI=merged_truncat_CI.sort_index()
    merged_truncat_CI.to_csv("merged_Truncated_CI_"+name+".csv", index=True)
    
os.popen("mv merged_Truncated_CI_* ../summary")
os.chdir("../summary")

### confirm that the rows and columns are matched
if sum(merged_gt_CI.columns == merged_truncat_CI.columns) != np.shape(merged_truncat_CI)[1]:
    sys.exit("Columns do not match, please double check the data")

if sum(merged_est_CI.index == merged_truncat_CI.index) != np.shape(merged_truncat_CI)[0]:
    sys.exit("Rows do not match, please double check the data")

### go to plot
# CI change box plot
def label_median(df, ax):
    medians=df.median()
    vertical_offset = max(medians.median(), 0.005)
    for xtick in ax.get_xticks():
        ax.text(xtick, medians[xtick] + vertical_offset, round(medians[xtick], 3), horizontalalignment='center',
                size='x-small', color='navy', weight='semibold')

# for each metagenome data, the col and row were already matched so we can substract/divide the dataframe element-wise directly
for line in query_list: #single plot for each metagenome
	name=line.split("/")[-1]
	print("Generating plot for "+name)
	temp_gt=pd.read_csv('merged_Ground_truth_CI_'+name+'.csv', header=0, index_col=0)
	temp_est=pd.read_csv('merged_Estimated_CI_'+name+'.csv', header=0, index_col=0)
	temp_trc=pd.read_csv('merged_Truncated_CI_'+name+'.csv', header=0, index_col=0)
	# boxplot
	fig1, axs = plt.subplots(3)
	sb.boxplot(data=temp_gt, color="red", ax=axs[0])
	sb.boxplot(data=temp_est, color="cyan", ax=axs[1])
	sb.boxplot(data=temp_trc, color="black", ax=axs[2])
	axs[0].set(ylabel="GroundTruth")
	axs[1].set(ylabel="Estimated")
	axs[2].set(ylabel="Truncated")
	fig1.tight_layout(rect=[0, 0.03, 1, 0.9])
	fig1.suptitle("Boxplot of 3 CIs")
	fig1.savefig("Boxplot_of_"+name+".png", dpi=200)
	plt.close(fig1)
	# swarmplot
	fig2, axs2 = plt.subplots(3)
	fig2.suptitle("Swarmplot of 3 CIs")
	sb.swarmplot(data=temp_gt, color="red", ax=axs2[0], size=3, linewidth=0.3)
	sb.swarmplot(data=temp_est, color="cyan", ax=axs2[1], size=3, linewidth=0.3)
	sb.swarmplot(data=temp_trc, color="black", ax=axs2[2], size=3, linewidth=0.3)
	axs2[0].set(ylabel="GroundTruth")
	axs2[1].set(ylabel="Estimated")
	axs2[2].set(ylabel="Truncated")
	fig2.tight_layout(rect=[0, 0.03, 1, 0.9])
	fig2.savefig("Swarmplot_of_"+name+".png", dpi=200)
	plt.close(fig2)
	### put 3 swarm together
	fig3, axs3 = plt.subplots()
	sb.swarmplot(data=temp_gt, color="red", size=3, linewidth=0.3, ax=axs3)
	sb.swarmplot(data=temp_est, color="cyan", size=3, linewidth=0.3, ax=axs3)
	sb.swarmplot(data=temp_trc, color="black", size=3, linewidth=0.3, ax=axs3)
	axs3.set(ylabel="CMash Index")
	legend_element = [Line2D([0],[0], marker='o', color='g', label='GroundTruth', markerfacecolor='red', markersize=5),
	                  Line2D([0],[0], marker='o', color='g', label='Estimated', markerfacecolor='cyan', markersize=5),
	                  Line2D([0],[0], marker='o', color='g', label='Truncated', markerfacecolor='black', markersize=5)]
	axs3.legend(handles=legend_element, loc='best')
	fig3.suptitle("Merged Swarmplot")
	fig3.savefig("Merged_swarmplot_of_"+name+".png", dpi=200)
	plt.close(fig3)
        # dif plot, temp_gt.min = 0.00003
	dif_trc_gt = temp_trc - temp_gt
	pdif_trc_gt = dif_trc_gt / (temp_gt+0.0000001)
	dif_est_gt = temp_est - temp_gt
	pdif_est_gt = dif_est_gt / (temp_gt+0.0000001)
	dif_trc_est = temp_trc - temp_est
	pdif_trc_est = dif_trc_est / (temp_est + 0.0000001)
	fig4, axs4 = plt.subplots(3,2, figsize=(16,12))
	sb.swarmplot(data=dif_trc_gt, ax=axs4[0,0], size=3, linewidth=0.3)
	sb.swarmplot(data=pdif_trc_gt, ax=axs4[0,1], size=3, linewidth=0.3)
	sb.swarmplot(data=dif_est_gt, ax=axs4[1,0], size=3, linewidth=0.3)
	sb.swarmplot(data=pdif_est_gt, ax=axs4[1,1], size=3, linewidth=0.3)
	sb.swarmplot(data=dif_trc_est, ax=axs4[2,0], size=3, linewidth=0.3)
	sb.swarmplot(data=pdif_trc_est, ax=axs4[2,1], size=3, linewidth=0.3)
	axs4[0,0].set(ylabel="Truncated vs GroundTruth")
	axs4[1,0].set(ylabel="Estimated vs GroundTruth")
	axs4[2,0].set(ylabel="Truncated vs Estimated")
	axs4[0,0].title.set_text("Change by number")
	axs4[0,1].title.set_text("Change by ratio")
	fig4.suptitle("Pairwise comparison")
	fig4.savefig("Pairwise_comparison_of_"+name+".png", dpi=200)
	plt.close(fig4)
	# boxplot of CI change: truncated vs groundtruth
	fig5, axs = plt.subplots(3, figsize=(9, 9))
	fig5.suptitle("Boxplot of CI dif")
	sb.boxplot(data=dif_trc_gt, color="red", ax=axs[0])
	label_median(dif_trc_gt, axs[0])
	sb.boxplot(data=dif_est_gt, color="cyan", ax=axs[1])
	label_median(dif_est_gt, axs[1])
	sb.boxplot(data=dif_trc_est, color="orange", ax=axs[2])
	label_median(dif_trc_est, axs[2])
	axs[0].set(ylabel="Truncated-Groundtruth")
	axs[1].set(ylabel="Estimated-Groundtruth")
	axs[2].set(ylabel="Truncated-Estimated")
	fig5.tight_layout(rect=[0, 0.03, 1, 0.9])
	fig5.savefig("Boxplot_of_CI_dif_"+name+".png", dpi=300)
	plt.close(fig5)







