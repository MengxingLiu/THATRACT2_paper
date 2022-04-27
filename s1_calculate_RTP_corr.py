import os
os.environ['FOR_DISABLE_CONSOLE_CTRL_HANDLER'] = '1'
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
import random, matplotlib, itertools, glob, platform, getpass
import scipy.stats
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
df = df.append(pd.read_csv(raw_csv_dir / "RTP_Profile_AL_07_fix.csv"))
con_fa_ana = pd.DataFrame(columns = ["subID", "TCK", "corr", "btw"])
for ana in ["AL_07", "AL_07_fix"]:
    tract_to_cal = df[df["analysis"]==ana].TCK.unique()
    df_ana_x = df[ (df["ses"] == "T01") & (df["analysis"] == ana)]
    df_ana_y = df[ (df["ses"] == "T02") & (df["analysis"] == ana)]
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
                                            "btw":"T01vsT02",
                                            "analysis":ana}, 
                                            ignore_index=True)
        

con_fa_ana = con_fa_ana.replace({"TCK":tract_dic} )
con_fa_ana = con_fa_ana.rename(columns={"subID":"SUBID"})
con_fa_ana["SUBID"] = con_fa_ana["SUBID"].fillna('') + con_fa_ana["SUBID2"].fillna('')
con_fa_ana.to_csv(git_dir / "correlation_fa_AL_07_withfix.csv", index=False)

correlation = pd.read_csv(git_dir / "correlation_fa.csv")

correlation.groupby("TCK").apply(np.mean)
correlation.groupby("TCK").describe().to_csv("correlation_description.csv")





df = pd.read_csv(raw_csv_dir / "RTP_Profile_compute.csv")
analysis = df.analysis.unique()
analysis = itertools.combinations(analysis,2)

## calculate RTP profile correlation

between = ["computation", "test-retest"]


inds = ['ad', 'cl', 'md', 'volume', 'curvature',
       'rd', 'fa', 'torsion']

con_fa_ana = pd.DataFrame()
for ana in analysis:
    for BTW in between:
        print(ana, " start")
        if BTW=="computation":
            
            df_ana_x = df[ (df["ses"] == "T01") & (df["analysis"] == ana[0])]
            df_ana_y = df[ (df["ses"] == "T01") & (df["analysis"] == ana[1])]
            btw = f"{ana[0]}vs{ana[1]}"
            df_ana = df_ana_x.merge(df_ana_y, how="inner",
                    on=["subID", "ind", "TCK", "ses"])
        elif ((BTW=="test-retest" ) & (ana == ("AL_07", "AL_08"))):

            df_ana_x = df[ (df["ses"] == "T01") & (df["analysis"] == ana[0])]
            df_ana_y = df[ (df["ses"] == "T02") & (df["analysis"] == ana[0])] 
            df_ana = df_ana_x.merge(df_ana_y, how="inner",
                    on=["subID", "ind", "TCK", "analysis"])
            btw = "test-retest_AL_07"
        else: 
            continue       
        tmp_df = pd.DataFrame()
        for i in inds:
            c = df_ana.groupby(["subID", "TCK"])[i+"_x", i+"_y"].corr().iloc[0::2,-1]
            c = c.reset_index()[["subID", "TCK", i+"_y"]]
            if i == inds[0]:
                tmp_df = c.copy()
            else: tmp_df[i+"_y"] = c[i+"_y"]
        tmp_df["btw"] = btw
        print(ana, " finish")
        con_fa_ana = con_fa_ana.append(tmp_df, ignore_index=True)

con_fa_ana_bk = con_fa_ana.copy()
con_fa_ana_bk.to_csv(raw_csv_dir / "con_fa_ana_raw.csv", index=False)
con_fa_ana = pd.read_csv(raw_csv_dir / "con_fa_ana_raw.csv")
con_filter = con_fa_ana[con_fa_ana["fa_y"]<0.7]
for row in con_filter.itertuples(index=True, name='Pandas'):
    sub = row.subID
    tck = row.TCK
    btw = row.btw
    fa = row.fa_y
    
    if "test-retest" in btw:
        df_x = df[((df.TCK==tck) & (df.subID==sub) & (df.ses=="T01") &
                    (df.analysis=="AL_07"))]
        df_y = df[((df.TCK==tck) & (df.subID==sub) & (df.ses=="T02") &
                    (df.analysis=="AL_07"))]
    else:
        df_x = df[((df.TCK==tck) & (df.subID==sub) & (df.ses=="T01") &
                    (df.analysis==btw.split("vs")[0]))]
        df_y = df[((df.TCK==tck) & (df.subID==sub) & (df.ses=="T01") &
                    (df.analysis==btw.split("vs")[1]))]
    b = scipy.stats.pearsonr(df_x["fa"], df_y["fa"][::-1])[0]
    if b < fa:
        
        continue
    else:
        print(sub, tck, btw, fa, b, "changing")
        for i in inds:
            b = scipy.stats.pearsonr(df_x[i], df_y[i][::-1])[0]
            con_fa_ana.loc[((con_fa_ana["TCK"]==tck) & 
                            (con_fa_ana.subID ==sub) &
                            (con_fa_ana.btw ==btw)), i+"_y"] = b

con_fa_ana = con_fa_ana.replace({"TCK":tract_dic} )
con_fa_ana = con_fa_ana.rename(columns={"subID":"SUBID"})
con_fa_ana.to_csv(git_dir / "correlation_all_index.csv", index=False)
