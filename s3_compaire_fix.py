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
raw_csv_dir = Path(f"{git_dir}/raw_csv")
# read tractparams to get the target label dictionary
tractparams = pd.read_csv(git_dir / "tractparams_AL_final_both_hemi.csv")
tract_dic = dict(zip(tractparams["slabel"], tractparams["roi2"]))


pairwise_TRT = pd.read_csv(raw_csv_dir / "pairwise_agreement.csv")
pairwise_TRT["btw"] = "T01vsT02"
pairwise_TRT["analysis"] = "AL_07"
pairwise_TRT_fix = pd.read_csv(raw_csv_dir / "pairwise_agreement_AL_07_fix.csv")
pairwise_TRT_fix["analysis"] = "AL_07_fix"
pairwise_TRT_fix["btw"] = "T01vsT02"
df = pd.concat([pairwise_TRT, pairwise_TRT_fix])
df = df.replace({"tract":tract_dic})
# comparie before and after fix
tract_to_plot = df[df["analysis"]=="AL_07_fix"].tract.unique()
fig, axes = plt.subplots()
sns.stripplot(x = "tract", y = "dice_voxels", 
                data = df[df["tract"].isin(tract_to_plot)],
                order = tract_to_plot, 
                hue = "analysis", ax = axes, alpha=0.5)
sns.pointplot(x = "tract", y = "dice_voxels", 
                order = tract_to_plot,
                data = df[df["tract"].isin(tract_to_plot)],
                hue = "analysis", alpha = 0.55, ax = axes, join=False)
plt.show()

# compare profile correlation
df = pd.read_csv(git_dir / "correlation_fa_AL_07_withfix.csv")
tract_to_plot = df[df["analysis"]=="AL_07_fix"].TCK.unique()
fig, axes = plt.subplots()
sns.stripplot(x = "TCK", y = "corr", 
                data = df[df["TCK"].isin(tract_to_plot)],
                order = tract_to_plot, 
                hue = "analysis", ax = axes, alpha=0.5)
sns.pointplot(x = "TCK", y = "corr", 
                order = tract_to_plot,
                data = df[df["TCK"].isin(tract_to_plot)],
                hue = "analysis", alpha = 0.55, ax = axes, join=False)
plt.show()



pairwise_com = pd.read_csv(raw_csv_dir / 
                        "pairwise_agreement_comAL_06vscomAL_07.csv")

pairwise_com_fix = pd.read_csv(raw_csv_dir / 
                        "pairwise_agreement_comAL_06_fixvscomAL_07_fix.csv")

df = pd.concat([pairwise_com, pairwise_com_fix])
df = df.replace({"tract":tract_dic})
# comparie before and after fix
tract_to_plot = df[df["btw"]=="comAL_06_fixvscomAL_07_fix"].tract.unique()
fig, axes = plt.subplots()
sns.stripplot(x = "tract", y = "dice_voxels", 
                data = df[df["tract"].isin(tract_to_plot)],
                order = tract_to_plot, 
                hue = "btw", ax = axes, alpha=0.5)
sns.pointplot(x = "tract", y = "dice_voxels", 
                order = tract_to_plot,
                data = df[df["tract"].isin(tract_to_plot)],
                hue = "btw", alpha = 0.55, ax = axes, join=False)
plt.show()


fig, axes = plt.subplots()
sns.stripplot(x = "tract", y = "density_correlation", 
                data = df[df["tract"].isin(tract_to_plot)],
                order = tract_to_plot, 
                hue = "btw", ax = axes, alpha=0.5)
sns.pointplot(x = "tract", y = "density_correlation", 
                order = tract_to_plot,
                data = df[df["tract"].isin(tract_to_plot)],
                hue = "btw", alpha = 0.55, ax = axes, join=False)
plt.show()


# profile correlation btw computations
df = pd.read_csv(git_dir / 
                    "correlation_fa_compute_withfix.csv")

tract_to_plot = df[df["btw"]=="AL_06_fixvsAL_07_fix"].TCK.unique()
fig, axes = plt.subplots()
sns.stripplot(x = "TCK", y = "corr", 
                data = df[df["TCK"].isin(tract_to_plot)],
                order = tract_to_plot, 
                hue = "btw", ax = axes, alpha=0.5)
sns.pointplot(x = "TCK", y = "corr", 
                order = tract_to_plot,
                data = df[df["TCK"].isin(tract_to_plot)],
                hue = "btw", alpha = 0.55, ax = axes, join=False)
plt.show()


pairwise_com = pd.read_csv(raw_csv_dir / 
                        "pairwise_agreement_comAL_06vscomAL_07.csv")

pairwise_com_fix = pd.read_csv(raw_csv_dir / 
                        "pairwise_agreement_comAL_06_fixvscomAL_07_fix.csv")


df[(df.TCK=="R_Area_25") &
     (df.btw=="AL_06_fixvsAL_07_fix")
     ].sort_values(by="corr")[:20]
df_pairwise = pd.concat([pairwise_com, pairwise_com_fix])
df_pairwise = df.replace({"tract":tract_dic})
df_pairwise[(df_pairwise.tract=="R_Area_25") & 
                (df_pairwise.btw=="comAL_06_fixvscomAL_07_fix")
                ].sort_values(by="dice_voxels")[:20]