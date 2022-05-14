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
    ax[0].set(title=f"{x} vs {y}")
    x="dice_voxels_T01vsT02_ln"
    sns.regplot(x=x,y=y, data=b, ax = ax[1], label = b[x].corr(b[y]))
    x="dice_voxels_comAL_all_average_ln"
    sns.regplot(x=x,y=y, data=b, ax = ax[1], label = b[x].corr(b[y]))
    ax[1].set(title=f"{x} vs {y}")
    x="density_correlation_T01vsT02_Fisher_z"
    sns.regplot(x=x,y=y, data=b, ax = ax[2], label = b[x].corr(b[y]))
    x="density_correlation_comAL_all_average_Fisher_z"
    sns.regplot(x=x,y=y, data=b, ax = ax[2], label = b[x].corr(b[y]))
    ax[2].set(title=f"{x} vs {y}")
    x="bundle_adjacency_voxels_T01vsT02"
    sns.regplot(x=x,y=y, data=b, ax = ax[3], label = b[x].corr(b[y]))
    x="bundle_adjacency_voxels_comAL_all_average"
    sns.regplot(x=x,y=y, data=b, ax = ax[3], label = b[x].corr(b[y]))
    ax[3].set(title=f"{x} vs {y}")
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

