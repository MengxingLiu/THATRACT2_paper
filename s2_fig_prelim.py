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


corr = pd.read_csv("correlation_fa.csv")
# plot Profile correlation between test-retest
corr_TRT = corr[corr["btw"]=="T01vsT02"]
corr_TRT_des = pd.read_csv("correlation_description.csv")

fig, axes = plt.subplots()
sns.stripplot(y="TCK", x="corr", order = corr_TRT_des.TCK.unique(), 
                data=corr_TRT, alpha = 0.5, ax = axes)
sns.pointplot(y="TCK", x="mean", order = corr_TRT_des.TCK.unique(), 
                data=corr_TRT_des, alpha = 0.55, ax = axes, join=False)
plt.xticks(rotation = -90)
plt.show()

# plot Profile correlation between AL06 AL07
corr_0607 = corr[corr["btw"]=="AL_06vsAL_07"]
corr_0607_des = corr_0607.groupby(["TCK"]).mean()
corr_0607_des = corr_0607_des.sort_values(by="corr")
corr_0607_des["TCK"] = corr_0607_des.index
fig, axes = plt.subplots()
sns.stripplot(y="TCK", x="corr", order = corr_0607_des.TCK.unique(), 
                data=corr_0607, alpha = 0.5, ax = axes)
sns.pointplot(y="TCK", x="corr", order = corr_0607_des.TCK.unique(), 
                data=corr_0607_des, alpha = 0.55, ax = axes, join=False)
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
sns.stripplot(y="TCK", x = "bundle_adjacency_voxels", data = pairwise_TRT, 
                order=order )
sns.pointplot(y="TCK", x = "bundle_adjacency_voxels", data = mean, 
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
sns.stripplot(y="TCK", x = "bundle_adjacency_voxels", data = pairwise_0607, 
                order=order )
sns.pointplot(y="TCK", x = "bundle_adjacency_voxels", data = mean, 
                order = order, join = False)
plt.xticks(rotation = -90)
#plt.tight_layout()
plt.show()


mean = mean.sort_values(by = "dice_voxels")
mean["TCK"]=mean.index
order = mean.TCK.unique()
fig, axes = plt.subplots()
sns.stripplot(y="TCK", x = "dice_voxels", data = pairwise_0607, 
                order=order )
sns.pointplot(y="TCK", x = "dice_voxels", data = mean, 
                order = order, join = False)
plt.xticks(rotation = -90)
plt.tight_layout()
plt.show()



# check the lowest values of specific tracts
corr = pd.read_csv("correlation_fa.csv")
pairwise = pd.read_csv(git_dir / "pairwise.csv")
pairwise_TRT = pairwise[pairwise["btw"]=="T01vsT02"]
# plot Profile correlation between test-retest
corr_TRT = corr[corr["btw"]=="T01vsT02"]
corr_TRT_des = pd.read_csv("correlation_description.csv")
# AL 12   50-130
corr_TRT[corr_TRT["TCK"] == "L_Area_25"].sort_values(by="corr")
pairwise_TRT[pairwise_TRT["TCK"] == "L_Area_25"][["SUBID","dice_voxels"]].sort_values(
                        by="dice_voxels")
# AL 5 nothing can be done, RTP profile problem   
corr_TRT[corr_TRT["TCK"] == "Left-MammillaryBody"].sort_values(by="corr")
pairwise_TRT[pairwise_TRT["TCK"] == "Left-MammillaryBody"][["SUBID","dice_voxels"]].sort_values(
                        by="dice_voxels")
# AL 25
corr_TRT[corr_TRT["TCK"] == "R_posterior_OFC_Complex"].sort_values(by="corr")
pairwise_TRT[pairwise_TRT["TCK"] == "R_posterior_OFC_Complex"][["SUBID","dice_voxels"]].sort_values(
                        by="dice_voxels")
# AL 3
corr_TRT[corr_TRT["TCK"] == "L_Entorhinal_Cortex"].sort_values(by="corr")
pairwise_TRT[pairwise_TRT["TCK"] == "L_Entorhinal_Cortex"][["SUBID","dice_voxels"]].sort_values(
                        by="dice_voxels")


# AL 23
corr_TRT[corr_TRT["TCK"] == "L_Area_13l"].sort_values(by="corr")
pairwise_TRT[pairwise_TRT["TCK"] == "L_Area_13l"][["SUBID","dice_voxels"]].sort_values(
                        by="dice_voxels")
# AL 14
corr_TRT[corr_TRT["TCK"] == "L_Area_s32"].sort_values(by="corr")
pairwise_TRT[pairwise_TRT["TCK"] == "L_Area_s32"][["SUBID","dice_voxels"]].sort_values(
                        by="dice_voxels")

# AL 24
corr_TRT[corr_TRT["TCK"] == "L_Orbital_Frontal_Complex"].sort_values(by="corr")
pairwise_TRT[pairwise_TRT["TCK"] == "L_Orbital_Frontal_Complex"][["SUBID","dice_voxels"]].sort_values(
                        by="dice_voxels")

# AL 3
corr_TRT[corr_TRT["TCK"] == "L_Entorhinal_Cortex"].sort_values(by="corr")
pairwise_TRT[pairwise_TRT["TCK"] == "L_Entorhinal_Cortex"][["SUBID","dice_voxels"]].sort_values(
                        by="dice_voxels")