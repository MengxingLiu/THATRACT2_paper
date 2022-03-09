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
tract_dic = dict(zip(tractparams["label"], tractparams["roi2"]))

# calculate results for test vs retest
df = pd.read_csv(raw_csv_dir / "RTP_Profile_AL_07.csv")

con_fa_ana = pd.DataFrame(columns = ["subID", "TCK", "corr", "btw"])

df_ana_x = df[ (df["ses"] == "T01") & (df["analysis"] == "AL_07")]
df_ana_y = df[ (df["ses"] == "T02") & (df["analysis"] == "AL_07")]
TCKS = df_ana_y.TCK.unique()
ses = df_ana_y.ses.unique()
SUBS = df_ana_y.subID.unique()
ind = 'fa'
# correlation of fa
for tck, sub in itertools.product(TCKS, SUBS):
    a = df_ana_x.loc[(df_ana_x["TCK"]==tck) & (df_ana_x["ses"]=="T01") 
                & (df_ana_x["subID"]==sub), "fa"]
    if len(a)==0 or a.isnull().values.any():
        continue
    b = df_ana_y.loc[(df_ana_y["TCK"]==tck) & (df_ana_y["ses"]=="T02") 
                & (df_ana_y["subID"]==sub), "fa"]
    if len(b)==0 or b.isnull().values.any():
        continue
    
    c, _ = scipy.stats.pearsonr(a,b)
    d, _ = scipy.stats.pearsonr(a[::-1],b)
    c = max(c,d)    
    print(tck, sub, "T01 vs T02")
    con_fa_ana = con_fa_ana.append({"subID":sub, "TCK":tck, 
                                          "corr":c, 
                                          "btw":"T01vsT02"}, 
                                         ignore_index=True)
    

con_fa_ana = con_fa_ana.replace({"TCK":tract_dic} )
con_fa_ana = con_fa_ana.rename(columns={"subID":"SUBID"})
con_fa_ana.to_csv(git_dir / "correlation_fa.csv", index=False)


correlation = pd.read_csv(git_dir / "correlation_fa.csv")
correlation.groupby("TCK").apply(np.mean)
correlation.groupby("TCK").describe().to_csv("correlation_description.csv")
