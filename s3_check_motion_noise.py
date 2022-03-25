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
head_motion = head_motion[[6,7,8]]
head_motion.columns = ["fwd", "SUBID", "ses"]
a = head_motion.groupby(["SUBID", "ses"]).mean().reset_index()
head_motion = a.pivot(index="SUBID", columns = "ses", values = "fwd").reset_index()
# average head motion across volumes

a.pivot(index=[7], columns=[8], values=[6]).reset_index()

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

scipy.stats.pearsonr(b["mean"], b["noise_T02"])
# -0.77 -0.79

corr_noise = pd.merge(corr_TRT, noise_TRT, on=["TCK", "SUBID"])
corr_noise_motion = pd.merge(corr_noise, head_motion, on="SUBID")
#check subject-wise correlation of specific tract
for i in corr_noise.TCK.unique():
    corr_noise_tck = corr_noise_motion[corr_noise_motion["TCK"]==i]
    r, p = scipy.stats.pearsonr(corr_noise_tck["T01_x"], corr_noise_tck["T01_y"])
    if r < -0.3 or r > 0.3:
        print(i,r, corr_noise_tck["corr"].mean(), "noise with motion")
    
    r, p = scipy.stats.pearsonr(corr_noise_tck["corr"], corr_noise_tck["T01_y"])
    if r < -0.3:
        print(i,r,p, corr_noise_tck["corr"].mean(), "fa correlation with motion")
        # only 8 tracts have <-0.3 correlation bewteen correlation and motion
        # not related to the lowest correlation TCKs

for i in corr_noise.SUBID.unique():
    corr_noise_sub = corr_noise[corr_noise["SUBID"]==i]
    r, p = scipy.stats.pearsonr(corr_noise_sub["corr"], corr_noise_sub["T01"])
    if r < -0.3:
        print(i,r, corr_noise_sub["corr"].mean())




# plot 
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

tmp = pd.merge(pairwise_TRT, noise_TRT, on = ["SUBID", "TCK"])
pairwise_noise_motion = pd.merge(tmp, head_motion, on=["SUBID"])


for i in pairwise_noise_motion.TCK.unique():
    pairwise_noise_motion_tck = pairwise_noise_motion[pairwise_noise_motion["TCK"]==i]
    # r, p = scipy.stats.pearsonr(corr_noise_tck["corr"], corr_noise_tck["T01_x"])
    # if r < -0.3:
    #     print(i,r, corr_noise_tck["corr"].mean())
    
    r, p = scipy.stats.pearsonr(pairwise_noise_motion_tck["dice_voxels"],
                         pairwise_noise_motion_tck["T01_y"])
    if r < -0.3:
        print(i,r,p, pairwise_noise_motion_tck["dice_voxels"].mean(), "dice with motion")
        # 10 tracts has <-0.3 correlation, not related to lower dice.



## test with fa correlation computational
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
sns.pointplot(y="TCK", x="T01", order = order, 
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
corr_noise_motion_com = pd.merge(corr_noise_compute, head_motion, on = "SUBID")

#check subject-wise correlation of specific tract
for i in corr_noise.TCK.unique():
    corr_noise_tck = corr_noise_motion_com[corr_noise_motion_com["TCK"]==i]
    r, p = scipy.stats.pearsonr(corr_noise_tck["T01_x"], corr_noise_tck["T01_y"])
    if r < -0.3 or r > 0.3:
        print(i,r, corr_noise_tck["corr"].mean(), "noise with motion")
    
    r, p = scipy.stats.pearsonr(corr_noise_tck["corr"], corr_noise_tck["T01_y"])
    if r < -0.3 or r > 0.3:
        print(i,r,p, corr_noise_tck["corr"].mean(), "fa correlation with motion")
        # only one tract has <-0.3 correlation of fa-motion


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

pairwise_noise_com = pd.merge(pairwise_0607, noise_pivot, on=["SUBID", "TCK"])
pairwise_noise_motion_com = pd.merge(pairwise_noise_com, head_motion, on="SUBID")

for i in pairwise_noise_com.TCK.unique():
    pairwise_noise_motion_tck = (
        pairwise_noise_motion_com[pairwise_noise_motion_com["TCK"]==i])
    # r, p = scipy.stats.pearsonr(corr_noise_tck["corr"], corr_noise_tck["T01_x"])
    # if r < -0.3:
    #     print(i,r, corr_noise_tck["corr"].mean())
    
    r = pairwise_noise_motion_tck["dice_voxels"].corr(
                 pairwise_noise_motion_tck["T01_y"])
    if r < -0.3 or r > 0.3:
        print(i,r, pairwise_noise_motion_tck["dice_voxels"].mean(), "dice with motion")
        # 1 tracts has <-0.3 correlation



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



pw_fa_ns_mt_st = pd.read_csv(git_dir / "pw_fa_ns_mt_st.csv")

# check correlation of test-retest correlation with noise
a = pw_fa_ns_mt_st.groupby("TCK").mean()
a = a.sort_values(by="fa_corr_06vs07")


plot_repro_factors(pw_fa_ns_mt_st,"fa_corr_T01vsT02", "noise_T01")
# -0.7866186472388043 -0.3682364176125297 if remove 10 highes and lowest
plot_repro_factors(pw_fa_ns_mt_st,"fa_corr_T01vsT02", "noise_T02")
# -0.7759863708526683 -0.46180884931222727 if remove 10 highes and lowest
plot_repro_factors(pw_fa_ns_mt_st,"fa_corr_T01vsT02", "noise_diff_T01-T02")
# 0.1180326456152848 0.23898675048453927 if remove 10 highes and lowest

plot_repro_factors(pw_fa_ns_mt_st,"tck_mean_T01", "noise_T01")
# -0.28565356632495487, 0.12796428929406803, if remove 10 highest and lowest 
# -0.11112394009215262 if remove mammilothalamic tract
plot_repro_factors(pw_fa_ns_mt_st,"tck_count_T01", "noise_T01")
# -0.08647096457729374, -0.29107895820172003, if remove 10 highest and lowest 
# -0.2797163980008855 if remove mammilothalamic tract


plot_repro_factors(pw_fa_ns_mt_st,"fa_corr_T01vsT02", "tck_mean_T01")
# 0.1290932172492908 -0.08636114611309718 if remove 10 highes and lowest
plot_repro_factors(pw_fa_ns_mt_st,"fa_corr_T01vsT02", "tck_count_T01")
# 0.3676500034555438 0.6934717535405337 if remove 10 highes and lowest


plot_repro_factors(pw_fa_ns_mt_st,"dice_voxels_T01vsT02", "noise_T01")
# -0.46863224868639336 -0.42998713825677454 if remove 10 highes and lowest
plot_repro_factors(pw_fa_ns_mt_st,"dice_voxels_T01vsT02", "noise_T02")
# -0.49079305367184406 -0.5138359427954071 if remove 10 highes and lowest
plot_repro_factors(pw_fa_ns_mt_st,"dice_voxels_T01vsT02", "noise_diff_T01-T02")
# 0.05715084780651699 0.37298867042959916 if remove 10 highes and lowest


plot_repro_factors(pw_fa_ns_mt_st,"dice_voxels_T01vsT02", "tck_mean_T01")
# -0.26310572058065407 0.16092644767706032 if remove 10 highes and lowest
plot_repro_factors(pw_fa_ns_mt_st,"dice_voxels_T01vsT02", "tck_count_T01")
# 0.7465002458634376 0.668889203349984 if remove 10 highes and lowest


plot_repro_factors(pw_fa_ns_mt_st,"fa_corr_06vs07", "noise_T01")
# -0.12397250548153714 -0.5728846646937854 if remove 10 highes and lowest


plot_repro_factors(pw_fa_ns_mt_st,"fa_corr_06vs07", "tck_mean_T01")
# -0.7401572126510043 0.5975722163072998 if remove 10 highes and lowest
plot_repro_factors(pw_fa_ns_mt_st,"fa_corr_06vs07", "tck_count_T01")
# 0.509771379182233 0.27489574508645587 if remove 10 highes and lowest


plot_repro_factors(pw_fa_ns_mt_st,"dice_voxels_comAL_06vscomAL_07", "noise_T01")
# -0.06672911036060773 -0.2337883790402507 if remove 10 highes and lowest

plot_repro_factors(pw_fa_ns_mt_st,"dice_voxels_comAL_06vscomAL_07", "tck_mean_T01")
# -0.787288825501465 0.5428474424217867 if remove 10 highes and lowest
plot_repro_factors(pw_fa_ns_mt_st,"dice_voxels_comAL_06vscomAL_07", "tck_count_T01")
# 0.6848783685323427 0.7756594415064574 if remove 10 highes and lowest


a["fa_corr_06vs07"].corr(a["noise_T01"])   # -0.12
#a[0:10]["fa_corr_06vs07"].corr(a[0:10]["noise_T01"])  # -0.54

a["fa_corr_06vs07"].corr(a["tck_mean_T01"])   # -0.74
a["fa_corr_06vs07"].corr(a["tck_count_T01"])   # 0.5



a["dice_voxels_comAL_06vscomAL_07"].corr(a["noise_T01"])   # -0.07
a[0:10]["dice_voxels_comAL_06vscomAL_07"].corr(a[0:10]["noise_T01"])  # -0.50
a["dice_voxels_comAL_06vscomAL_07"].corr(a["tck_mean_T01"])   # -0.78
a["dice_voxels_comAL_06vscomAL_07"].corr(a["tck_count_T01"])   # 0.68
a[5:-5]["dice_voxels_comAL_06vscomAL_07"].corr(a[5:-5]["tck_mean_T01"])

plot_repro_factors(pw_fa_ns_mt_st,"dice_voxels_comAL_06vscomAL_07", "noise_T01")
sns.set_style("darkgrid")
def plot_repro_factors(data, repro, factor):
    a = data.groupby("TCK").mean()
    a = a.sort_values(by=repro).reset_index()
    r = a[repro].corr(a[factor]) 
    r_10 = a[10:-10][repro].corr(a[10:-10][factor])
    b = a[~a["TCK"].str.contains("Mamm")]
    r_no_mam = b[repro].corr(b[factor])
    print(f"{r}, {r_10}, if remove 10 highest and lowest ")
    print(f"{r_no_mam} if remove mammilothalamic tract")
    fig, axes = plt.subplots(figsize=(38,20))
    axes2 = axes.twiny()
    # fa correlation with noise
    sns.stripplot(y="TCK", x=repro, order = a["TCK"], 
                    data=data, alpha = 0.5, ax = axes)
    sns.pointplot(y="TCK", x=repro, order = a["TCK"], 
                    data=a, alpha = 0.55, ax = axes, 
                    join=False, color="Blue")
    sns.pointplot(y="TCK", x=factor, order = a["TCK"], 
                    data=a, alpha = 0.55, ax = axes2, 
                    join=False, color="Green")
    axes.set_title(f"correlation between {repro} and {factor}/ r = {r}")
    plt.xticks(rotation = -90)
    #plt.show() #plt.show()








a["fa_corr_T01vsT02"].corr(a["noise_T01"])   # -0.79
a = a.sort_values(by="fa_corr_T01vsT02")
a[0:10]["fa_corr_T01vsT02"].corr(a[0:10]["noise_T01"])  # -0.87

# subject-wise correlation
for i in pw_fa_ns_mt_st.TCK.unique():
    corr_noise_tck = pw_fa_ns_mt_st[pw_fa_ns_mt_st["TCK"]==i]
    r = corr_noise_tck["noise_T01"].corr(corr_noise_tck["motion_T01"])
    if r < -0.3 or r > 0.3:
        print(i,r, corr_noise_tck["fa_corr_06vs07"].mean(), "noise with motion")
    
    r = corr_noise_tck["fa_corr_06vs07"].corr(corr_noise_tck["motion_T01"])
    if r < -0.3 or r > 0.3:
        print(i,r, corr_noise_tck["fa_corr_06vs07"].mean(), "fa correlation with motion")
    
    r = corr_noise_tck["fa_corr_06vs07"].corr(corr_noise_tck["noise_T01"])
    if r < -0.3 or r > 0.3:
        print(i,r, corr_noise_tck["fa_corr_06vs07"].mean(), "fa correlation with noise")
    
    # conclusion: noise has no relation with motion
    # motion has subject-wise correlation with fa correlation only one tract, 8B_lateral
    # noist has positive relation with fa correlation on 6 tracts


# tract-wise correlation
for i in pw_fa_ns_mt_st.SUBID.unique():
    corr_noise_sub = pw_fa_ns_mt_st[pw_fa_ns_mt_st["SUBID"]==i]
    r = corr_noise_sub["noise_T01"].corr(corr_noise_sub["motion_T01"])
    if r < -0.3 or r > 0.3:
        print(i,r, corr_noise_sub["fa_corr_06vs07"].mean(), "noise with motion")
    
    r = corr_noise_sub["fa_corr_06vs07"].corr(corr_noise_sub["motion_T01"])
    if r < -0.3 or r > 0.3:
        print(i,r, corr_noise_sub["fa_corr_06vs07"].mean(), "fa correlation with motion")
    
    r = corr_noise_sub["fa_corr_06vs07"].corr(corr_noise_sub["noise_T01"])
    if r < -0.3 or r > 0.3:
        print(i,r, corr_noise_sub["fa_corr_06vs07"].mean(), "fa correlation with noise")
    # conclusion: noise has no relation with motion
    # noise has subject-wise correlation with fa correlation only one tract, ToTPilot69
    


# subject-wise correlation
for i in pw_fa_ns_mt_st.TCK.unique():
    corr_noise_tck = pw_fa_ns_mt_st[pw_fa_ns_mt_st["TCK"]==i]
    r = corr_noise_tck["fa_corr_T01vsT02"].corr(corr_noise_tck["noise_T01"])
    if r < -0.3 or r > 0.3:
        print(i,r, corr_noise_tck["fa_corr_T01vsT02"].mean(), "fa correlation with noise")
        # many tracts has <-0.3 correlation between fa correlation with noise

    r = corr_noise_tck["dice_voxels_T01vsT02"].corr(corr_noise_tck["noise_T01"])
    if r < -0.3 or r > 0.3:
        print(i,r, corr_noise_tck["dice_voxels_T01vsT02"].mean(), "dice with noise")
        


# tract-wise correlation
for i in pw_fa_ns_mt_st.SUBID.unique():
    corr_noise_sub = pw_fa_ns_mt_st[pw_fa_ns_mt_st["SUBID"]==i]
    r = corr_noise_sub["fa_corr_T01vsT02"].corr(corr_noise_sub["noise_T01"])
    if r < -0.3 or r > 0.3:
        print(i,r, corr_noise_sub["fa_corr_T01vsT02"].mean(), "fa correlation with noise")
        # many tracts has <-0.3 correlation between fa correlation with noise

    r = corr_noise_sub["dice_voxels_T01vsT02"].corr(corr_noise_sub["noise_T01"])
    if r < -0.3 or r > 0.3:
        print(i,r, corr_noise_sub["dice_voxels_T01vsT02"].mean(), "dice with noise")
        