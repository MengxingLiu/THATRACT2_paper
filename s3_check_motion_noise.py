from telnetlib import AO
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

# load head motion, the column 6 is FWD
head_motion = pd.read_csv(raw_csv / "head_motion_all.csv", header=None)

# calculate FWD percentage > 1mm for each subject
motion_summary = pd.DataFrame()
thrld = 0.5; 
SUBS = head_motion[7].unique(); SES = head_motion[8].unique()
for sub, ses in itertools.product(SUBS, SES):
    fwd = head_motion[((head_motion[7]==sub) & 
                (head_motion[8]==ses))][6]
    if len(fwd) == 0: continue
    fwd_ex = fwd[fwd>thrld]
    fwd_percent = len(fwd_ex)/len(fwd)
    motion_summary = motion_summary.append(
                    {"sub": sub, "ses": ses, 
                    "thrld":thrld, "percent":fwd_percent}, ignore_index=True)
    if fwd_percent > 0.05:
        print(sub, ses, fwd_percent)


tract_to_check = ['L_Area_25', 'R_Area_25', 'L_posterior_OFC_Complex',
       'R_posterior_OFC_Complex', 'L_Area_13l', 'R_Area_13l',
       'L_Area_s32', 'R_Area_s32', 'L_Orbital_Frontal_Complex',
       'R_Orbital_Frontal_Complex', 'Left-MammillaryBody',
       'Right-MammillaryBody']


pairwise_TRT = pd.read_csv(raw_csv / "pairwise_agreement_TRT.csv")
fa_corr = pd.read_csv(git_dir / "correlation_fa.csv")
fa_corr_TRT = fa_corr[fa_corr["btw"]=="T01vsT02"]

# check visually in excel
pairwise_TRT[pairwise_TRT["TCK"].isin(tract_to_check)].to_csv(
                    git_dir / "check_head_motion.csv", index=False )
fa_corr_TRT[fa_corr_TRT.TCK.isin(tract_to_check)].to_csv(
                    git_dir / "tmp.csv", index=False)

motion_summary.to_csv(git_dir/ "tmp_0.5.csv", index = False)


# check noise with fa test-retest correlation

corr = pd.read_csv("correlation_fa.csv")
# plot Profile correlation between test-retest
corr_TRT = corr[corr["btw"]=="T01vsT02"]
corr_TRT_des = pd.read_csv("correlation_description.csv")
noise = pd.read_csv(git_dir / "noise_clean.csv")
noise_TRT = noise[noise.SUBID.isin(corr_TRT.SUBID.unique())]
noise_TRT = noise_TRT.pivot(index = ["SUBID", "TCK"], columns = "ses", values="noise")
noise_TRT = noise_TRT.reset_index()
noise_des = noise_TRT.groupby("TCK").mean().reset_index()


fig, axes = plt.subplots()
axes2 = axes.twiny()
sns.stripplot(y="TCK", x="corr", order = corr_TRT_des.TCK.unique(), 
                data=corr_TRT, alpha = 0.5, ax = axes)
sns.pointplot(y="TCK", x="mean", order = corr_TRT_des.TCK.unique(), 
                data=corr_TRT_des, alpha = 0.55, ax = axes, join=False)
sns.pointplot(y="TCK", x="T01", order = corr_TRT_des.TCK.unique(), 
                data=noise_des, alpha = 0.55, ax = axes2, 
                join=False, color="Green")

sns.pointplot(y="TCK", x="T02", order = corr_TRT_des.TCK.unique(), 
                data=noise_des, alpha = 0.55, ax = axes2, 
                join=False, color="Red")

plt.xticks(rotation = -90)
plt.show()

b = pd.merge(corr_TRT_des, noise_des, on="TCK")

scipy.stats.pearsonr(b["mean"], b["T02"])
corr_noise = pd.merge(corr_TRT, noise_TRT, on=["TCK", "SUBID"])
#check subject-wise correlation of specific tract
for i in corr_noise.TCK.unique():
    corr_noise_tck = corr_noise[corr_noise["TCK"]==i]
    r, p = scipy.stats.pearsonr(corr_noise_tck["corr"], corr_noise_tck["T02"])
    if r < -0.3:
        print(i,r, corr_noise_tck["corr"].mean())

for i in corr_noise.SUBID.unique():
    corr_noise_sub = corr_noise[corr_noise["SUBID"]==i]
    r, p = scipy.stats.pearsonr(corr_noise_sub["corr"], corr_noise_sub["T01"])
    if r < -0.3:
        print(i,r, corr_noise_sub["corr"].mean())


corr_a = corr_TRT[corr_TRT["TCK"]=="L_Area_25"].sort_values(by="SUBID")
noise_a = noise_TRT[noise_TRT["TCK"]=="L_Area_25"].sort_values(by="SUBID")
scipy.stats.pearsonr(corr_a["corr"], noise_a["T01"])
# check noise with dice test-retest
pairwise_TRT = pd.read_csv(raw_csv / "pairwise_agreement_TRT.csv")
mean = pairwise_TRT.groupby("TCK").mean()

mean = mean.sort_values(by = "dice_voxels")
mean = mean.reset_index()
order = mean.TCK.unique()
fig, axes = plt.subplots()
axes2 = axes.twiny()
sns.stripplot(y="TCK", x = "dice_voxels", data = pairwise_TRT, 
                order=order, ax=axes )
sns.pointplot(y="TCK", x = "dice_voxels", data = mean, 
                order = order, join = False, ax = axes)
sns.pointplot(y="TCK", x="T01", order = order, 
                data=noise_des, alpha = 0.55, ax = axes2, 
                join=False, color="Green")

sns.pointplot(y="TCK", x="T02", order = order, 
                data=noise_des, alpha = 0.55, ax = axes2, 
                join=False, color="Red")

plt.xticks(rotation = -90)
plt.tight_layout()
plt.show()
a = pd.merge(mean, noise_des, on="TCK")

scipy.stats.pearsonr(a["dice_voxels"], a["T02"])
-0.50, -0.49


## test with reproducibility fa correlation
corr_0607 = corr[corr["btw"]=="AL_06vsAL_07"]
corr_0607_des = corr_0607.groupby(["TCK"]).mean()
corr_0607_des = corr_0607_des.sort_values(by="corr")
corr_0607_des = corr_0607_des.reset_index()

noise_pivot = noise.pivot(index = ["SUBID", "TCK"], columns = "ses", values="noise")
noise_pivot = noise_pivot.reset_index()
noise_des = noise_pivot.groupby("TCK").mean().reset_index()
order = corr_0607_des.TCK.unique()
fig, axes = plt.subplots()
axes2 = axes.twiny()

sns.stripplot(y="TCK", x="corr", order = order, 
                data=corr_0607, alpha = 0.5, ax = axes)
sns.pointplot(y="TCK", x="corr", order = order, 
                data=corr_0607_des, alpha = 0.55, ax = axes, join=False)
sns.pointplot(y="TCK", x="noise", order = order, 
                data=noise_des, alpha = 0.55, ax = axes2, 
                join=False, color="Green")

# sns.pointplot(y="TCK", x="T02", order = order, 
#                 data=noise_des, alpha = 0.55, ax = axes2, 
#                 join=False, color="Red")
plt.xticks(rotation = -90)
plt.show()

a = pd.merge(corr_0607_des, noise_des, on="TCK")

scipy.stats.pearsonr(a["corr"][0:10], a["noise"][0:10])
scipy.stats.pearsonr(a["corr"], a["noise"])
# try drop mammilory body
a = a.drop(a[a["TCK"].str.contains("MammillaryBody")].index)
-0.12

corr_noise_compute = pd.merge(corr_0607, noise_pivot,
                             on=["SUBID", "TCK"])

#check subject-wise correlation of specific tract
for i in corr_noise_compute.TCK.unique():
    corr_noise_compute_tck = corr_noise_compute[corr_noise_compute["TCK"]==i]
    r, p = scipy.stats.pearsonr(corr_noise_compute_tck["corr"], corr_noise_compute_tck["T01"])
    if r < -0.2:
        print(i,r, corr_noise_compute_tck["corr"].mean())

for i in corr_noise_compute.SUBID.unique():
    corr_noise_compute_sub = corr_noise_compute[corr_noise_compute["SUBID"]==i]
    r, p = scipy.stats.pearsonr(corr_noise_compute_sub["corr"], corr_noise_compute_sub["T01"])
    if r < -0.3:
        print(i,r, corr_noise_compute_sub["corr"].mean())



## test with reproducibility dice
pairwise = pd.read_csv(git_dir / "pairwise.csv")
pairwise_0607 = pairwise[pairwise["btw"]=="comAL_06vscomAL_07"]
mean = pairwise_0607.groupby("TCK").mean()
mean = mean.sort_values(by = "dice_voxels")
mean = mean.reset_index()
order = mean.TCK.unique()
fig, axes = plt.subplots()
axes2 = axes.twiny()
sns.stripplot(y="TCK", x = "dice_voxels", data = pairwise_0607, 
                order=order, ax=axes )
sns.pointplot(y="TCK", x = "dice_voxels", data = mean, 
                order = order, join = False, ax=axes)
sns.pointplot(y="TCK", x="noise", order = order, 
                data=noise_des, alpha = 0.55, ax = axes2, 
                join=False, color="Green")
        
plt.xticks(rotation = -90)
#plt.tight_layout()
plt.show()


a = pd.merge(mean, noise_des, on="TCK")

scipy.stats.pearsonr(a["dice_voxels"], a["noise"])
scipy.stats.pearsonr(a["dice_voxels"][0:15], a["noise"][0:15])
a[["TCK","dice_voxels", "noise"]][0:15]
-0.07, -0.07



for i in corr_0607.TCK.unique():
    c = corr_0607[corr_0607["TCK"]==i]
    noise_com =noise[noise.SUBID.isin(c.SUBID.unique())]
    noise_com = noise_com[noise_com["TCK"]==i]
    noise_com = noise_com[noise_com["ses"]=="T01"]
    c = c.sort_values(by="SUBID")
    noise_com = noise_com.sort_values(by="SUBID")
    r, p = scipy.stats.pearsonr(c["corr"], noise_com.noise)
    if r < -0.3:
        print(i,r )