# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 14:05:49 2021

@author: lmengxing
"""
import pandas as pd
import seaborn as sns
import itertools
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy
import numpy as np
import random
import matplotlib.gridspec as gridspec
""" calculate results for compute vs recompute
tractDic = {"LKN27":"L_OR_05", "LKN28":"R_OR_05", "LKN29":"L_OR_1", "LKN30":"R_OR_1",
            "LKN31":"L_AR_belt-3", "LKN32":"L_AR_belt-4",
            "LKN33":"L_AR_A1-3", "LKN34": "L_AR_A1-4",
            "LKN35":"R_AR_belt-3", "LKN36":"R_AR_belt-4",
            "LKN37":"R_AR_A1-3", "LKN38":"R_AR_A1-4",
            "LKN39":"L_DT-3", "LKN40":"L_DT-4",
            "LKN41":"R_DT-3", "LKN42":"R_DT-4",
            "LKN43":"L_MR_M1", "LKN44":"R_MR_M1",
            "LKN45":"L_MR_dlPreM", "LKN46":"R_MR_dlPreM", 
            "LKN47":"L_MR_S1", "LKN48":"R_MR_S1"      }


list1 = [f"{i:02d}" for i in range(1,11)]
list1 = range(1,11)
a = itertools.combinations(list1, 2)
df = pd.read_csv("RTP_Profile_compute.csv")

con_fa_ana = pd.DataFrame(columns = ["subID", "TCK", "corr", "btw"])

for i in a:
    
    df_ana_x = df[ (df["ses"] == "T01") & (df["analysis"] == i[0])]
    df_ana_y = df[ (df["ses"] == "T01") & (df["analysis"] == i[1])]
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
        b = df_ana_y.loc[(df_ana_y["TCK"]==tck) & (df_ana_y["ses"]=="T01") 
                    & (df_ana_y["subID"]==sub), "fa"]
        if len(b)==0 or b.isnull().values.any():
            continue
        
        c, _ = scipy.stats.pearsonr(a,b)
        d, _ = scipy.stats.pearsonr(a[::-1],b)
        c = max(c,d)    
        print(tck, sub, f"{i[0]} and {i[1]}")
        con_fa_ana = con_fa_ana.append({"subID":sub, "TCK":tck, 
                                              "corr":c, 
                                              "btw":f"{i[0]}vs{i[1]}"}, 
                                             ignore_index=True)
        
con_fa_ana["TCK"] = con_fa_ana.TCK.map(lambda x : "L"+x)
con_fa_ana = con_fa_ana.replace({"TCK":tractDic} )
con_fa_ana = con_fa_ana.replace({"TCK":newlabel})
con_fa_ana.to_csv("correlation_fa_compute.csv", index=False)

" calculate descriptions"

correlation = pd.read_csv("correlation_fa_compute.csv")
pairwise = pd.read_csv("pairwise_agreement_compute.csv")

pairwise = pairwise.replace({"tract":tractDic} )

tck_to_plot = ["L_OR_05", "R_OR_05", "L_AR_A1-4", "R_AR_A1-4", 
               "L_MR_M1", "R_MR_M1", "L_DT-4", "R_DT-4"]
profile = profile[profile["TCK"].isin(tck_to_plot)]
correlation = correlation[correlation["TCK"].isin(tck_to_plot)]
correlation["TCK"] = correlation["TCK"].map(lambda x : x[0:4])
                                                 
pairwise = pairwise[pairwise["tract"].isin(tck_to_plot)]
pairwise["tract"] = pairwise["tract"].map(lambda x : x[0:4])
pairwise["bundle_adjacency_voxels"] = pairwise["bundle_adjacency_voxels"].astype(float)
pairwise["dice_voxels"] = pairwise["dice_voxels"].astype(float)
pairwise["density_correlation"] = pairwise["density_correlation"].astype(float)

correlation.groupby("TCK").apply(np.mean)
correlation.groupby("TCK").describe().to_csv("correlation_compute_description.csv")
pairwise.groupby("tract").describe().to_csv("pairwise_agreement_compute_description.csv")
         
"""

""" calculate results for test vs retest
csv_dir = "F:\TESTDATA\THATRACT_paper\csv"

df = pd.read_csv(f"{csv_dir}\RTP_Profile.csv")

con_fa_ana = pd.DataFrame(columns = ["subID", "TCK", "corr", "btw"])

df_ana_x = df[ (df["ses"] == "T01") & (df["analysis"] == 1)]
df_ana_y = df[ (df["ses"] == "T02") & (df["analysis"] == 1)]
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
    
con_fa_ana["TCK"] = con_fa_ana.TCK.map(lambda x : "L"+x)
con_fa_ana = con_fa_ana.replace({"TCK":tractDic} )
con_fa_ana = con_fa_ana.replace({"TCK":newlabel})
con_fa_ana.to_csv("correlation_fa.csv", index=False)  

" calculate descriptions "

correlation = pd.read_csv("correlation_fa.csv")
pairwise = pd.read_csv("pairwise_agreement.csv")
pairwise = pairwise.replace({"tract":tractDic} )

tck_to_plot = ["L_OR_05", "R_OR_05", "L_AR_A1-4", "R_AR_A1-4", 
               "L_MT_M1", "R_MT_M1", "L_DT-4", "R_DT-4"]
profile = profile[profile["TCK"].isin(tck_to_plot)]
correlation = correlation[correlation["TCK"].isin(tck_to_plot)]
correlation["TCK"] = correlation["TCK"].map(lambda x : x[0:4])
                                                 
pairwise = pairwise[pairwise["tract"].isin(tck_to_plot)]
pairwise["tract"] = pairwise["tract"].map(lambda x : x[0:4])
correlation.groupby("TCK").apply(np.mean)

correlation.groupby("TCK").describe().to_csv("correlation_description.csv")
pairwise.groupby("tract").describe().to_csv("pairwise_agreement_description.csv")
           
"""

""" plot Fig. 2 
csv_dir = "F:\TESTDATA\GIT\THATRACT_paper\csv"
profile = pd.read_csv(f"{csv_dir}\RTP_Profile_compute.csv")
correlation = pd.read_csv(f"{csv_dir}\correlation_fa_compute.csv")
pairwise = pd.read_csv(f"{csv_dir}\pairwise_agreement_compute.csv")

# change tck labels
profile["TCK"] = profile.TCK.map(lambda x : "L"+x)
profile = profile.replace({"TCK":tractDic} )

pairwise = pairwise.replace({"tract":tractDic} )

tck_to_plot = ["L_OR_05", "R_OR_05", "L_AR_A1-4", "R_AR_A1-4", 
               "L_MR_M1", "R_MR_M1", "L_DT-4", "R_DT-4"]
profile = profile[profile["TCK"].isin(tck_to_plot)]

correlation = correlation.replace({"TCK":{"L_MT_M1":"L_MR_M1","R_MT_M1":"R_MR_M1"}})
correlation = correlation[correlation["TCK"].isin(tck_to_plot)]
correlation["TCK"] = correlation["TCK"].map(lambda x : x[0:4])
                                                 
pairwise = pairwise[pairwise["tract"].isin(tck_to_plot)]
pairwise["tract"] = pairwise["tract"].map(lambda x : x[0:4])

pairwise["bundle_adjacency_voxels"] = pairwise["bundle_adjacency_voxels"].astype(float)
pairwise["dice_voxels"] = pairwise["dice_voxels"].astype(float)
pairwise["density_correlation"] = pairwise["density_correlation"].astype(float)

sns.set_style("darkgrid")
fig2 = plt.figure(constrained_layout=True)
gs = fig2.add_gridspec(6, 2)
f2_ax1 = fig2.add_subplot(gs[0:3, 0])
palette = sns.color_palette("Paired")
order = ["L_OR","R_OR","L_AR", "R_AR", "L_MR", "R_MR", "L_DT", "R_DT"]
# FA profile 

tmp = profile[(profile["analysis"].isin([1,2])) & (profile["TCK"]=='L_MR_M1') &
              (profile["ses"]=="T01")]
tmp["analysis"] = tmp["analysis"].map(lambda x : f"compute_0{str(x)}")
f2_ax1 = sns.lineplot(data=tmp, x="ind", y="fa", hue="analysis",
                   ci="sd", style = "analysis", palette = ["Grey", "Green"] )
f2_ax1.set_title('gs[0, :]')
f2_ax2 = fig2.add_subplot(gs[3:6, 0])

# FA correlation dist
f2_ax2 = sns.stripplot(x="TCK", y = "corr", data = correlation,
                       order = order, palette = palette, rasterized=True)
f2_ax2.set(xlabel="fiber group label")

f2_ax3 = fig2.add_subplot(gs[0:2, 1])
f2_ax3 = sns.stripplot(x="tract", y = "bundle_adjacency_voxels", data = pairwise,
                    order= order, palette = palette,
                    rasterized=True)
f2_ax3.set(xticklabels=[])
f2_ax3.set(xlabel=None)

f2_ax4 = fig2.add_subplot(gs[2:4, 1])
f2_ax4 = sns.stripplot(x="tract", y = "dice_voxels", data = pairwise,
                       order = order, palette = palette, rasterized=True)
f2_ax4.set(xticklabels=[])
f2_ax4.set(xlabel=None)

f2_ax5 = fig2.add_subplot(gs[4:6, 1])
f2_ax5 = sns.stripplot(x="tract", y = "density_correlation", data = pairwise,
                       order = order, palette = palette, rasterized=True)
f2_ax5.set(xlabel="fiber group label")

fig_dir = r"F:\TESTDATA\GIT\THATRACT_paper\figures"
fig2.set_size_inches(9.88,4.93)
fig2.savefig(f"{fig_dir}\Fig2_computational_new.svg", dpi=300, bbox_inches='tight')
"""

 # plot Fig.3 
csv_dir = r"F:\TESTDATA\GIT\THATRACT_paper\csv"
profile = pd.read_csv(f"{csv_dir}\RTP_Profile.csv")
correlation = pd.read_csv(f"{csv_dir}\correlation_fa.csv")
pairwise = pd.read_csv(f"{csv_dir}\pairwise_agreement.csv")

# change tck labels
profile["TCK"] = profile.TCK.map(lambda x : "L"+x)
profile = profile.replace({"TCK":tractDic} )
pairwise = pairwise.replace({"tract":tractDic} )

tck_to_plot = ["L_OR_05", "R_OR_05", "L_AR_A1-4", "R_AR_A1-4", 
               "L_MR_M1", "R_MR_M1", "L_DT-4", "R_DT-4"]
profile = profile[profile["TCK"].isin(tck_to_plot)]


profile = profile[profile["subID"].isin(profile[profile["ses"]=="T02"].subID.unique())]
correlation = correlation.replace({"TCK":{"L_MT_M1":"L_MR_M1","R_MT_M1":"R_MR_M1"}})
correlation = correlation[correlation["TCK"].isin(tck_to_plot)]
correlation["TCK"] = correlation["TCK"].map(lambda x : x[0:4])
                                                 
pairwise = pairwise[pairwise["tract"].isin(tck_to_plot)]
pairwise["tract"] = pairwise["tract"].map(lambda x : x[0:4])

pairwise["bundle_adjacency_voxels"] = pairwise["bundle_adjacency_voxels"].astype(float)
pairwise["dice_voxels"] = pairwise["dice_voxels"].astype(float)
pairwise["density_correlation"] = pairwise["density_correlation"].astype(float)

fig3 = plt.figure(constrained_layout=True)
gs = fig3.add_gridspec(6, 2)
f3_ax1 = fig3.add_subplot(gs[0:3, 0])

# FA profile 
#tmp = profile[(profile["TCK"]=='L_MR_M1') & (profile["subID"]=="S038")]
tmp = profile[(profile["TCK"]=='L_MR_M1')]
tmp = tmp[tmp["subID"].isin(tmp[tmp["ses"]=="T02"].subID.unique())]
f3_ax1 = sns.lineplot(data=tmp, x="ind", y="fa", hue="ses",
                   style="ses", ci="sd", palette = ["Grey", "Green"] )
# FA correlation dist
f3_ax2 = fig3.add_subplot(gs[3:6, 0])
f3_ax2 = sns.stripplot(x="TCK", y = "corr", data = correlation,
                       order = order, palette = palette, rasterized=True)
f3_ax2.set(xlabel="fiber group label")

f3_ax3 = fig3.add_subplot(gs[0:2, 1])

f3_ax3 = sns.stripplot(x="tract", y = "bundle_adjacency_voxels",
                       data = pairwise, order = order, palette = palette,
                       rasterized=True)
f3_ax3.set(xticklabels=[])
f3_ax3.set(xlabel=None)

f3_ax4 = fig3.add_subplot(gs[2:4, 1])
f3_ax4 = sns.stripplot(x="tract", y = "dice_voxels", data = pairwise, 
                       order = order, palette = palette, rasterized=True)
f3_ax4.set(xticklabels=[])
f3_ax4.set(xlabel=None)

f3_ax5 = fig3.add_subplot(gs[4:6, 1])
f3_ax5 = sns.stripplot(x="tract", y = "density_correlation", data = pairwise,
                       order = order, palette = palette, rasterized=True)
f3_ax5.set(xlabel="fiber group label")

fig3.set_size_inches(9.88,4.93)
fig3.savefig(f"{fig_dir}\Fig3_test-retest_new.svg", dpi=300, format = "svg", bbox_inches='tight')
"""