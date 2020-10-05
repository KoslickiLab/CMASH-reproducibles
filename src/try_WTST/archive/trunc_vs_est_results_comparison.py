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
input_range="20-60-5"
os.chdir("./download")


start=input_range.split("-")[0]
end=input_range.split("-")[1]
gap=input_range.split("-")[2]



# generate a list of variance number / percentage for 2 matched files
def generate_list(file1, file2, symmatric=False):
    df1 = pd.read_csv(file1, header=0, index_col=0)
    df2 = pd.read_csv(file2, header=0, index_col=0) #used as base
    dif = df1 - df2
    dif_ratio = dif / (df2+0.00001)
    # if the matrix is self-matched (symmatric, diagnol=1), remove those unuseful records
    if symmatric:
        dif = dif.where(np.triu(np.ones(dif.shape), 1).astype(np.bool))
        dif_ratio = dif_ratio.where(np.triu(np.ones(dif_ratio.shape), 1).astype(np.bool))
        df1 = df1.where(np.triu(np.ones(df1.shape), 1).astype(np.bool))
        df2 = df2.where(np.triu(np.ones(df2.shape), 1).astype(np.bool))
    value_dif = dif.stack().to_numpy()
    ratio_dif = dif_ratio.stack().to_numpy()
    value_df1 = df1.stack().to_numpy()
    value_df2 = df2.stack().to_numpy()
    return value_df1, value_df2, value_dif, ratio_dif


# read csv and merge the results into dataframe for plotting
df_trunc = pd.DataFrame()
df_est = pd.DataFrame()
df_dif = pd.DataFrame()
df_dif_ratio = pd.DataFrame()
for k in range(int(start), int(end), int(gap)):
    temp_colname = 'k='+str(k)
    v1, v2, vdif, rdif = generate_list("trunc_WJI_k"+str(k)+".csv", "est_WJI_k"+str(k)+".csv", symmatric=True)
    df_trunc[temp_colname] = v1
    df_est[temp_colname] = v2
    df_dif[temp_colname] = vdif
    df_dif_ratio[temp_colname] = rdif



# generate plots
name="high_simi_data"
print("Generating plot for "+name)

# boxplot
fig1, axs = plt.subplots(3)
sb.boxplot(data=df_trunc, color="red", ax=axs[0])
sb.boxplot(data=df_est, color="cyan", ax=axs[1])
sb.boxplot(data=df_dif, color="black", ax=axs[2])
axs[0].set(ylabel="Trunc_WJI")
axs[1].set(ylabel="Est_WJI")
axs[2].set(ylabel="Trunc-Est")
fig1.tight_layout(rect=[0, 0.03, 1, 0.9])
fig1.suptitle("Boxplot of WJI")
fig1.savefig("Boxplot_of_" + name + ".png", dpi=200)
plt.close(fig1)

### merged swam plot
fig3, axs3 = plt.subplots()
sb.swarmplot(data=df_trunc, color="red", size=3, linewidth=0.3, ax=axs3)
sb.swarmplot(data=df_est, color="cyan", size=3, linewidth=0.3, ax=axs3)
#sb.swarmplot(data=temp_trc, color="black", size=3, linewidth=0.3, ax=axs3)
axs3.set(ylabel="WJI value")
legend_element = [Line2D([0], [0], marker='o', color='g', label='Trunc', markerfacecolor='red', markersize=5),
                  Line2D([0], [0], marker='o', color='g', label='Est', markerfacecolor='cyan', markersize=5),
                  Line2D([0], [0], marker='o', color='g', label='GT(not now)', markerfacecolor='black', markersize=5)]
axs3.legend(handles=legend_element, loc='best')
fig3.suptitle("Merged Swarmplot")
fig3.savefig("Merged_swarmplot_of_" + name + ".png", dpi=200)
plt.close(fig3)



# dif plot, temp_gt.min = 0.00003
fig4, axs4 = plt.subplots(2, figsize=(16,12))
sb.swarmplot(data=df_dif, ax=axs4[0], size=3, linewidth=0.3)
sb.swarmplot(data=df_dif_ratio, ax=axs4[1], size=3, linewidth=0.3)
axs4[0].set(ylabel="Truncated vs Est")
axs4[1].title.set_text("Change by ratio")
fig4.suptitle("Pairwise comparison")
fig4.savefig("Pairwise_comparison_of_"+name+".png", dpi=200)
plt.close(fig4)

def label_median(df, ax):
    medians=df.median()
    vertical_offset = max(medians.median(), 0.005)
    for xtick in ax.get_xticks():
        ax.text(xtick, medians[xtick] + vertical_offset, round(medians[xtick], 3), horizontalalignment='center',
                size='x-small', color='navy', weight='semibold')

# boxplot of CI change: truncated vs groundtruth
fig5, axs = plt.subplots(1, figsize=(9, 9))
fig5.suptitle("Boxplot of WJI dif")
sb.boxplot(data=df_dif, color="red")
label_median(df_dif, axs)
axs.set(ylabel="Truncated-Groundtruth")
fig5.tight_layout(rect=[0, 0.03, 1, 0.9])
fig5.savefig("Boxplot_of_WJI_dif_" + name + ".png", dpi=300)
plt.close(fig5)