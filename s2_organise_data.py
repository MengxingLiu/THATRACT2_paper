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


## load pairwise_TRT
pairwise_TRT = pd.read_csv(raw_csv_dir / "pairwise_agreement.csv")
pairwise_TRT["btw"] = "T01vsT02"
pairwise_TRT["analysis"] = "AL_07"

## load pairwise computational
pairwise_0607 = pd.read_csv(raw_csv_dir / "pairwise_agreement_comAL_06vscomAL_07.csv")
pairwise = pd.concat([pairwise_TRT, pairwise_0607])
pairwise = pairwise.replace({"tract":tract_dic})
pairwise = pairwise.rename(columns={"tract":"TCK"})
pairwise["bundle_adjacency_voxels"] = pairwise["bundle_adjacency_voxels"].astype(float)
pairwise["dice_voxels"] = pairwise["dice_voxels"].astype(float)
pairwise["density_correlation"] = pairwise["density_correlation"].astype(float)
pairwise.to_csv(git_dir / "pairwise.csv", index=False)

pairwise = pairwise[["bundle_adjacency_voxels", "dice_voxels", 
                'density_correlation', 'TCK', 'SUBID', 'btw']]

pairwise = pairwise.pivot(index=["SUBID", "TCK"], columns = "btw")
pairwise.columns = ["_".join(i) for i in pairwise.columns.to_flat_index()]
pairwise = pairwise.reset_index()

noise = pd.read_csv(raw_csv_dir / "noise.csv")
noise = noise.drop(noise[noise["SUBID"]=="SUBID"].index )
noise["noise"] = noise["noise"].apply(lambda x: x.lstrip("b'").rstrip(r" \n'"))
noise = noise.replace({"TCK":tract_dic})
noise = noise.pivot(index=["SUBID", "TCK"], columns = "ses", values = "noise")
noise = noise.reset_index()
noise = noise.rename(columns = {"T01": "noise_T01", "T02": "noise_T02"})
noise.to_csv(git_dir / "noise_clean.csv", index=False)

# concatenate pairwise and noise
pairwise_noise = pd.merge(pairwise, noise,  how = "outer", on = ["SUBID", "TCK"])

# load head motion, the column 6 is FWD
head_motion = pd.read_csv(raw_csv_dir / "head_motion_all.csv", header=None)
head_motion = head_motion[[6,7,8]]
head_motion.columns = ["fwd", "SUBID", "ses"]
# average head motion across volumes
a = head_motion.groupby(["SUBID", "ses"]).mean().reset_index()
head_motion = a.pivot(index="SUBID", columns = "ses", values = "fwd").reset_index()
head_motion.columns = ["SUBID", "motion_T01", "motion_T02"]

# concatenate pairwise, noise and head_motion
pw_ns_mt = pd.merge(pairwise_noise, head_motion, how = "outer", on = ["SUBID"])


# load fa correlation
fa_corr = pd.read_csv(git_dir / "correlation_fa.csv")
fa_corr = fa_corr.pivot(index=["SUBID", "TCK"], columns= "btw", values = "corr")
fa_corr = fa_corr.reset_index()
fa_corr.columns = ["SUBID", "TCK", "fa_corr_06vs07", "fa_corr_T01vsT02"]

# concatenate pairwise, noise, head_motion and fa correlation
pw_fa_ns_mt = pd.merge(pw_ns_mt, fa_corr, how = "outer", on = ["SUBID", "TCK"])

# load tck stats
tckstats = pd.read_csv(git_dir / "tckstats_AL_07.csv")
tckstats = tckstats.pivot(index = ["SUBID", "TCK"], columns = "ses")
tckstats.columns = ["_".join(i) for i in tckstats.columns.to_flat_index()]
tckstats = tckstats.reset_index()

# concatenate pairwise, noise, head_motion, fa correlation and tckstats
pw_fa_ns_mt_st = pd.merge(pw_fa_ns_mt, tckstats, 
                        how="outer", on = ["SUBID", "TCK"])

pw_fa_ns_mt_st.to_csv(git_dir / "pw_fa_ns_mt_st.csv", index=False)

"""

# plot pairwise agreement bewtwee test-retest
pairwise_TRT = pd.read_csv(raw_csv_dir / "pairwise_agreement.csv")
pairwise_TRT["btw"] = "T01vsT02"
pairwise_TRT["analysis"] = "AL_07"
pairwise_TRT_fix = pd.read_csv(raw_csv_dir / "pairwise_agreement_AL_07_fix.csv")
pairwise_TRT_fix["analysis"] = "AL_07_fix"
pairwise_TRT_fix["btw"] = "T01vsT02"
df = pd.concat([pairwise_TRT, pairwise_TRT_fix])
df = df.replace({"tract":tract_dic})
"""
# tckstats organize
tckstats = pd.read_csv(raw_csv_dir / 
                "tckstats_AL_07.csv")

tckstats["stats"] = tckstats.stats.map(lambda x:x.lstrip("b'").rstrip(r"\\n'"))
tckstats[["tck_mean", "tck_std", "tck_min", "tck_max", "tck_count"]] = (
            tckstats["stats"].str.split(' ',4,  expand=True)
)
tckstats = tckstats.replace({"TCK": tract_dic})
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

# calculate description
pairwise.groupby(["btw", "TCK"]).describe().to_csv("pairwise_description.csv")


# organize noise file
noise = pd.read_csv(raw_csv_dir / "noise.csv")
noise = noise.drop(noise[noise["SUBID"]=="SUBID"].index )
noise["noise"] = noise["noise"].apply(lambda x: x.lstrip("b'").rstrip(r" \n'"))
noise = noise.replace({"TCK":tract_dic})
noise.to_csv(git_dir / "noise_clean.csv", index=False)