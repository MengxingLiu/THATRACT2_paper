from telnetlib import AO
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy,random, matplotlib, itertools, glob, os, platform, getpass, re
import numpy as np
import matplotlib.gridspec as gridspec
from pathlib import Path

if getpass.getuser() == "mengxing":
    git_dir = Path("/home/mengxing/GIT/THATRACT2_paper")
elif getpass.getuser() == "lmengxing":
    if platform.system() == "Linux":
        git_dir = Path("/bcbl/home/home_g-m/lmengxing/TESTDATA/GIT/THATRACT2_paper")
    elif platform.system() == "Windows":
        git_dir = Path("F:\TESTDATA\GIT\THATRACT2_paper")
raw_csv = Path(f"{git_dir}/raw_csv")



pw_fa_ns_mt_st = pd.read_csv(raw_csv / "pw_fa_ns_mt_st.csv")
# drop those computational stats
for col in pw_fa_ns_mt_st.columns:
    if re.match(".*AL.*AL", col):
        del pw_fa_ns_mt_st[col]

# calculate group mean of each correlation pair
# first plot reproducibility with noise
data = pw_fa_ns_mt_st.copy()
a = data.groupby("TCK").mean()
a = a.reset_index()
b = a[~a["TCK"].str.contains("Mamm")]
factors = ["noise_T01", "tck_mean_T01", "tck_count_T01"]
for i in factors:
    y = i
    plt.close()
    sns.set_style("darkgrid")
    fig, ax = plt.subplots(2,2)
    ax = ax.flatten()
    x="test-retest_AL_07_fa_Fisher_Z"
    sns.regplot(x=x,y=y, data=b, ax = ax[0], label = b[x].corr(b[y]))
    x="AL_all_fa_Fisher_Z"; 
    sns.regplot(x=x,y=y, data=b, ax = ax[0], label = b[x].corr(b[y]))
    ax[0].set(title=f"FA correlation Fisher_Z vs {y}")
    x="dice_voxels_T01vsT02_ln"
    sns.regplot(x=x,y=y, data=b, ax = ax[1], label = b[x].corr(b[y]))
    x="dice_voxels_comAL_all_average_ln"
    sns.regplot(x=x,y=y, data=b, ax = ax[1], label = b[x].corr(b[y]))
    ax[1].set(title=f"ln Dice index vs {y}")
    x="density_correlation_T01vsT02_Fisher_z"
    sns.regplot(x=x,y=y, data=b, ax = ax[2], label = b[x].corr(b[y]))
    x="density_correlation_comAL_all_average_Fisher_z"
    sns.regplot(x=x,y=y, data=b, ax = ax[2], label = b[x].corr(b[y]))
    ax[2].set(title=f" density correlation Fisher_z vs {y}")
    x="bundle_adjacency_voxels_T01vsT02"
    sns.regplot(x=x,y=y, data=b, ax = ax[3], label = b[x].corr(b[y]))
    x="bundle_adjacency_voxels_comAL_all_average"
    sns.regplot(x=x,y=y, data=b, ax = ax[3], label = b[x].corr(b[y]))
    ax[3].set(title=f"bundle adjacency vs {y}")
    for a in ax:
        a.legend()
    # plt.legend()
    plt.show()


corr_pairs = {"tck_mean_T01": "noise_T01", "tck_count_T01": "noise_T01",
             "test-retest_AL_07_fa_Fisher_Z": "noise_T01",
             "test-retest_AL_07_fa_Fisher_Z": "noise_T02",
             "test-retest_AL_07_fa_Fisher_Z": "tck_mean_T01",
             "test-retest_AL_07_fa_Fisher_Z": "tck_count_T01",
             "dice_voxels_T01vsT02_ln": "noise_T01",
             "dice_voxels_T01vsT02_ln": "tck_mean_T01",
             "dice_voxels_T01vsT02_ln": "tck_count_T01",
             "density_correlation_T01vsT02_Fisher_z": "noise_T01",
             "density_correlation_T01vsT02_Fisher_z": "tck_mean_T01",
             "density_correlation_T01vsT02_Fisher_z": "tck_count_T01",
             "bundle_adjacency_voxels_T01vsT02": "noise_T01",
             "bundle_adjacency_voxels_T01vsT02": "tck_mean_T01",
             "bundle_adjacency_voxels_T01vsT02": "tck_count_T01",
             "AL_all_fa_Fisher_Z": "noise_T01",
             "AL_all_fa_Fisher_Z": "tck_mean_T01",
             "AL_all_fa_Fisher_Z": "tck_count_T01",
             "dice_voxels_comAL_all_average_ln": "noise_T01",
             "dice_voxels_comAL_all_average_ln": "tck_mean_T01",
             "dice_voxels_comAL_all_average_ln": "tck_count_T01",
             "density_correlation_comAL_all_average_Fisher_z": "noise_T01",
             "density_correlation_comAL_all_average_Fisher_z": "tck_mean_T01",
             "density_correlation_comAL_all_average_Fisher_z": "tck_count_T01",                           
             "bundle_adjacency_voxels_comAL_all_average": "noise_T01",
             "bundle_adjacency_voxels_comAL_all_average": "tck_mean_T01",
             "bundle_adjacency_voxels_comAL_all_average": "tck_count_T01"
               }


repro = ["test-retest_AL_07_fa_Fisher_Z", "dice_voxels_T01vsT02_ln", 
        "density_correlation_T01vsT02_Fisher_z", "bundle_adjacency_voxels_T01vsT02", 
        "AL_all_fa_Fisher_Z", "dice_voxels_comAL_all_average_ln", 
        "density_correlation_comAL_all_average_Fisher_z", 
        "bundle_adjacency_voxels_comAL_all_average"]
corr_heat = pd.DataFrame()
for key, value in itertools.product(repro, factors):
    r = b[key].corr(b[value])
    corr_heat = corr_heat.append({"repro":key, value:r}, ignore_index=True)
corr_heat = corr_heat.groupby("repro").max()

TRT = ['test-retest_AL_07_fa_Fisher_Z',
        'bundle_adjacency_voxels_T01vsT02',
       'density_correlation_T01vsT02_Fisher_z',
       'dice_voxels_T01vsT02_ln']
COM = ['AL_all_fa_Fisher_Z', 'bundle_adjacency_voxels_comAL_all_average',
       'density_correlation_comAL_all_average_Fisher_z',
       'dice_voxels_comAL_all_average_ln']
fig, ax = plt.subplots(1,2)
cmap = "RdBu_r"
cmap = sns.diverging_palette(0, 230, 90, 60, as_cmap=True)
xlabels = ["noise", "length", "number"]
ylabels = ["FA correlation", "adjacency", "density", "Dice"]
sns.heatmap(corr_heat.loc[TRT], annot=True, fmt=".2f", cmap=cmap, 
                    center=0, vmin=-1, vmax=1,
                    xticklabels = xlabels,yticklabels=ylabels,
                    ax = ax[0])

sns.heatmap(corr_heat.loc[COM], annot=True, fmt=".2f", cmap=cmap, 
                    center=0, vmin=-1, vmax=1,
                    xticklabels = xlabels,yticklabels=ylabels,
                    ax = ax[1])
plt.show()

plt.close()
sns.set(style='white', font_scale=1.6)
fig, axes = plt.subplots(4,3)
for ind, ax in zip(itertools.product(TRT, factors), axes.flatten()):
    corr_r = corr_heat.loc[ind[0]][ind[1]]
    corr_text = f"{corr_r:2.2f}".replace("0.", ".")
    if 0.209 < abs(corr_r) < 0.27:
        sig_text = "*"
    elif 0.27 < abs(corr_r) < 0.34:
        sig_text = "**"
    elif 0.24 < abs(corr_r) < 1: 
        sig_text= "***"
    else: sig_text = ""
    marker_size = abs(corr_r)*10000
    ax.set_axis_off()
    ax.scatter([.5],[.5], marker_size, [corr_r],
            alpha = 0.6, vmin=-1,vmax=1, cmap = cmap)
    font_size = abs(corr_r) * 40 + 5       
    ax.annotate(corr_text, [.5, .5,],  xycoords="axes fraction",
                ha='center', va='center', fontsize=font_size)
    ax.annotate(sig_text, [.5, .7,],  xycoords="axes fraction",
                ha='center', va='center',
                color='red', fontsize=20)
plt.show()

sns.set_style("darkgrid")
fig2 = plt.figure(constrained_layout=True)
gs = fig2.add_gridspec(6, 2)
f2_ax1 = fig2.add_subplot(gs[0:3, 0])
palette = sns.color_palette("Paired")
order = ["L_OR","R_OR","L_AR", "R_AR", "L_MR", "R_MR", "L_DT", "R_DT"]
# FA profile 

tmp = profile[(profile["analysis"].isin([1,2])) & (profile["TCK"]=='L_MR_M1') &
              (profile["ses"]=="T01")]
tmp["analysis"] = tmp["analysis"].map(lambda x : f"compute_0{str(x)}")
f2_ax1 = sns.lineplot(data=tmp, x="ind", y="fa", hue="analysis",
                   ci="sd", style = "analysis", palette = ["Grey", "Green"] )
f2_ax1.set_title('gs[0, :]')
f2_ax2 = fig2.add_subplot(gs[3:6, 0])

# FA correlation dist
f2_ax2 = sns.stripplot(x="TCK", y = "corr", data = correlation,
                       order = order, palette = palette, rasterized=True)
f2_ax2.set(xlabel="fiber group label")

f2_ax3 = fig2.add_subplot(gs[0:2, 1])
f2_ax3 = sns.stripplot(x="tract", y = "bundle_adjacency_voxels", data = pairwise,
                    order= order, palette = palette,
                    rasterized=True)
f2_ax3.set(xticklabels=[])
f2_ax3.set(xlabel=None)

f2_ax4 = fig2.add_subplot(gs[2:4, 1])
f2_ax4 = sns.stripplot(x="tract", y = "dice_voxels", data = pairwise,
                       order = order, palette = palette, rasterized=True)
f2_ax4.set(xticklabels=[])
f2_ax4.set(xlabel=None)

f2_ax5 = fig2.add_subplot(gs[4:6, 1])
f2_ax5 = sns.stripplot(x="tract", y = "density_correlation", data = pairwise,
                       order = order, palette = palette, rasterized=True)
f2_ax5.set(xlabel="fiber group label")

fig_dir = r"F:\TESTDATA\GIT\THATRACT_paper\figures"
fig2.set_size_inches(9.88,4.93)
fig2.savefig(f"{fig_dir}\Fig2_computational_new.svg", dpi=300, bbox_inches='tight')