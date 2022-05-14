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




pw_fa_ns_mt_st = pd.read_csv(raw_csv / "pw_fa_ns_mt_st.csv")
for col in pw_fa_ns_mt_st.columns:
    if re.match(".*AL.*AL", col):
        del pw_fa_ns_mt_st[col]
# check correlation of test-retest correlation with noise
a = pw_fa_ns_mt_st.groupby("TCK").mean()
a = a.sort_values(by="test-retest_AL_07_fa")

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
rep = ["test-retest_AL_07_fa_Fisher_Z", "dice_voxels_T01vsT02_ln", 
        "density_correlation_T01vsT02_Fisher_z","bundle_adjacency_voxels_T01vsT02", 
        "AL_all_fa_Fisher_Z","dice_voxels_comAL_all_average_ln",
         "density_correlation_comAL_all_average_Fisher_z",
         "bundle_adjacency_voxels_comAL_all_average"]
factors = ["noise_T01", "noise_T02","tck_mean_T01" ,"tck_count_T01"]
for key, value in itertools.product(rep,factors):
    print(key, value)
    plot_repro_factors(pw_fa_ns_mt_st, key, value)

""" old commands
plot_repro_factors(pw_fa_ns_mt_st,"test-retest_AL_07_fa", "noise_T01")
# -0.7866186472388043 -0.3682364176125297 if remove 10 highes and lowest
plot_repro_factors(pw_fa_ns_mt_st,"test-retest_AL_07_fa_Fisher_Z", "noise_T01")
# -0.7719620804445716, -0.5420327708800353, if remove 10 highest and lowest
# -0.7322329086835533 if remove mammilothalamic tract

plot_repro_factors(pw_fa_ns_mt_st,"test-retest_AL_07_fa", "noise_T02")
# -0.7759863708526683 -0.46180884931222727 if remove 10 highes and lowest
plot_repro_factors(pw_fa_ns_mt_st,"test-retest_AL_07_fa_Fisher_Z", "noise_T02")
# -0.7963196858608802, -0.6245193851144515, if remove 10 highest and lowest 
# -0.7690002074294058 if remove mammilothalamic tract


plot_repro_factors(pw_fa_ns_mt_st,"test-retest_AL_07_fa", "noise_diff_T01-T02")
# 0.1180326456152848 0.23898675048453927 if remove 10 highes and lowest

plot_repro_factors(pw_fa_ns_mt_st,"tck_mean_T01", "noise_T01")
# -0.28565356632495487, 0.12796428929406803, if remove 10 highest and lowest 
# -0.11112394009215262 if remove mammilothalamic tract
plot_repro_factors(pw_fa_ns_mt_st,"tck_count_T01", "noise_T01")
# -0.08647096457729374, -0.29107895820172003, if remove 10 highest and lowest 
# -0.2797163980008855 if remove mammilothalamic tract


plot_repro_factors(pw_fa_ns_mt_st,"test-retest_AL_07_fa", "tck_mean_T01")
# 0.1290932172492908 -0.08636114611309718 if remove 10 highes and lowest
plot_repro_factors(pw_fa_ns_mt_st,"test-retest_AL_07_fa_Fisher_Z", "tck_mean_T01")
# 0.03875062208253553, 0.006052888355612519, if remove 10 highest and lowest 
# -0.0961565510908329 if remove mammilothalamic tract
plot_repro_factors(pw_fa_ns_mt_st,"test-retest_AL_07_fa", "tck_count_T02")
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

"""

plot_repro_factors(pw_fa_ns_mt_st,"dice_voxels_comAL_06vscomAL_07", "noise_T01")
sns.set_style("darkgrid")
def plot_repro_factors(data, repro, factor):
    plt.close()
    data = data[((data[repro].notna()) & (data[factor].notna()))]
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
        
