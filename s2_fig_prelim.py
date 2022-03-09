import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy,random, matplotlib, itertools, glob, os, platform, getpass
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
sns.set_style("darkgrid")
# plot Profile correlation between test-retest
corr_TRT = pd.read_csv("correlation_fa.csv")
corr_TRT_des = pd.read_csv("correlation_description.csv")
fig, axes = plt.subplots()
sns.stripplot(y="TCK", x="corr", order = corr_TRT_des.TCK.unique(), 
                data=corr_TRT, alpha = 0.5, ax = axes)
sns.pointplot(y="TCK", x="mean", order = corr_TRT_des.TCK.unique(), 
                data=corr_TRT_des, alpha = 0.55, ax = axes, join=False)
plt.xticks(rotation = -90)
plt.show()

### plot pairwise agreement bewtwee test-retest
# bundle adjacency 
pairwise = pd.read_csv(git_dir / "pairwise.csv")
pairwise_TRT = pairwise[pairwise["btw"]=="T01vsT02"]
mean = pairwise_TRT.groupby("TCK").mean()
mean = mean.sort_values(by = "bundle_adjacency_voxels")
mean["TCK"]=mean.index
order = mean.TCK.unique()
fig, axes = plt.subplots()
sns.stripplot(x="TCK", y = "bundle_adjacency_voxels", data = pairwise_TRT, 
                order=order )
sns.pointplot(x="TCK", y = "bundle_adjacency_voxels", data = mean, 
                order = order, join = False)
plt.xticks(rotation = -90)
plt.tight_layout()
plt.show()

# dice 
mean = mean.sort_values(by = "dice_voxels")
mean["TCK"]=mean.index
order = mean.TCK.unique()
fig, axes = plt.subplots()
sns.stripplot(y="TCK", x = "dice_voxels", data = pairwise_TRT, 
                order=order )
sns.pointplot(y="TCK", x = "dice_voxels", data = mean, 
                order = order, join = False)
plt.xticks(rotation = -90)
plt.tight_layout()
plt.show()



### plot pairwise agreement bewtwee computational
# bundle adjacency 
pairwise = pd.read_csv(git_dir / "pairwise.csv")
pairwise_0607 = pairwise[pairwise["btw"]=="comAL_06vscomAL_07"]
mean = pairwise_0607.groupby("TCK").mean()
mean = mean.sort_values(by = "bundle_adjacency_voxels")
mean["TCK"]=mean.index
order = mean.TCK.unique()
fig, axes = plt.subplots()
sns.stripplot(x="TCK", y = "bundle_adjacency_voxels", data = pairwise_0607, 
                order=order )
sns.pointplot(x="TCK", y = "bundle_adjacency_voxels", data = mean, 
                order = order, join = False)
plt.xticks(rotation = -90)
#plt.tight_layout()
plt.show()
