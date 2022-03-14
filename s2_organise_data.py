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

# plot pairwise agreement bewtwee test-retest
pairwise_TRT = pd.read_csv(raw_csv_dir / "pairwise_agreement.csv")
pairwise_TRT["btw"] = "T01vsT02"
pairwise_TRT["analysis"] = "AL_07"
pairwise_TRT_fix = pd.read_csv(raw_csv_dir / "pairwise_agreement_AL_07_fix.csv")
pairwise_TRT_fix["analysis"] = "AL_07_fix"
pairwise_TRT_fix["btw"] = "T01vsT02"
df = pd.concat([pairwise_TRT, pairwise_TRT_fix])
df = df.replace({"tract":tract_dic})

# tckstats organize
tckstats = pd.read_csv(raw_csv_dir / 
                "tckstats_AL_07.csv")

tckstats["stats"] = tckstats.stats.map(lambda x:x.lstrip("b'").rstrip(r"\\n'"))
tckstats[["mean", "std", "min", "max", "count"]] = (
            tckstats["stats"].str.split(' ',4,  expand=True)
)
tckstats.to_csv(git_dir / 
            "tckstats_AL_07.csv", index=False)
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


pairwise_0607 = pd.read_csv(raw_csv_dir / "pairwise_agreement_comAL_06vscomAL_07.csv")
pairwise = pd.concat([pairwise_TRT, pairwise_0607])
pairwise = pairwise.replace({"tract":tract_dic})
pairwise = pairwise.rename(columns={"tract":"TCK"})
pairwise["bundle_adjacency_voxels"] = pairwise["bundle_adjacency_voxels"].astype(float)
pairwise["dice_voxels"] = pairwise["dice_voxels"].astype(float)
pairwise["density_correlation"] = pairwise["density_correlation"].astype(float)
pairwise.to_csv(git_dir / "pairwise.csv", index=False)

# calculate description
pairwise.groupby(["btw", "TCK"]).describe().to_csv("pairwise_description.csv")

