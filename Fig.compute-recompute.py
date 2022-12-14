# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 14:05:49 2021

@author: lmengxing
"""
import os
os.environ['FOR_DISABLE_CONSOLE_CTRL_HANDLER'] = '1'
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy,random, matplotlib, itertools, glob, os, platform, getpass
import numpy as np
import matplotlib.gridspec as gridspec
from pathlib import Path
import matplotlib.ticker as ticker
plt.rcParams['svg.fonttype'] = 'none'
if getpass.getuser() == "mengxing":
    git_dir = Path("/home/mengxing/GIT/THATRACT2_paper")
elif getpass.getuser() == "lmengxing":
    if platform.system() == "Linux":
        git_dir = Path("/bcbl/home/home_g-m/lmengxing/TESTDATA/GIT/THATRACT2_paper")
    elif platform.system() == "Windows":
        git_dir = Path("F:\TESTDATA\GIT\THATRACT2_paper")
raw_csv = Path(f"{git_dir}/raw_csv")
fig_dir = Path(f"{git_dir}/figures")


profile = pd.read_csv( raw_csv / "RTP_Profile_AL_compute.csv")
correlation = pd.read_csv( raw_csv / "correlation_all_index.csv")
pairwise = pd.read_csv( raw_csv / "pairwise_agreement_all_final.csv")
pairwise = pairwise[(~(pairwise.dice_voxels==0) & 
                    ~(pairwise.density_correlation==0))]
# select computational results
#correlation = correlation[correlation["btw"]=="AL_06vsAL_07"]
#pairwise = pairwise[pairwise.btw=="comAL_06vscomAL_07"]



def plot_fig(btw, profile=profile, correlation=correlation, 
                pairwise=pairwise, tck_to_plot=tck_to_plot, pallete=palette):
    
    profile = profile[profile["TCK"].isin(tck_to_plot)]
    correlation = correlation[correlation["TCK"].isin(tck_to_plot)]
    pairwise = pairwise[pairwise["TCK"].isin(tck_to_plot)]
    size=7
    if btw == "compute":
        tmp = profile[(profile["analysis"].isin(["AL_06","AL_07"])) & 
                        (profile["TCK"]==tck_to_plot[0]) &
                        (profile["ses"]=="T01")]
        #tmp["analysis"] = tmp["analysis"].map(lambda x : f"compute_{x}")
        hue_pro = "analysis"
        correlation = correlation[  ~(correlation["btw"]=="test-retest_AL_07")]
        pairwise = pairwise[  ~(pairwise["btw"]=="T01vsT02")]
        labels = ["1st computation", "2nd computation"]
    elif btw == "test-retest":
        tmp = profile[(profile["TCK"]==tck_to_plot[0])]
        tmp = tmp[tmp["subID"].isin(tmp[tmp["ses"]=="T02"].subID.unique())]
        tmp = tmp[tmp["analysis"]=="AL_07"]
        hue_pro = "ses"
        labels = ["test", "retest"]
        correlation = correlation[  (correlation["btw"]=="test-retest_AL_07")]
        pairwise = pairwise[(pairwise["btw"]=="T01vsT02")]
    
    order = list(tck_to_plot).copy()
    xtick = [x[2:] for x in order]
    new_palette = palette.copy() * (round(len(order)/len(palette))+1)
    i = 2
    
    while i < len(order):
        # add space between bundle groups
        order.insert(i, "")
        order.insert(i+1, "")
        new_palette.insert(i, palette[0]);new_palette.insert(i, palette[0])
        i += (3+1)

    # FA profile 
    print(order)
    dodge = 0.2

    sns.set_style("darkgrid")
    sns.set(font_scale = 2)
    fig2 = plt.figure(constrained_layout=True)
    gs = fig2.add_gridspec(6, 2)
    f2_ax1 = fig2.add_subplot(gs[0:3, 0])
    f2_ax1 = sns.lineplot(data=tmp, x="ind", y="fa", hue=hue_pro,
                    ci="sd", style = hue_pro, palette = ["Grey", "Green"] )
    f2_ax1.set(ylabel="FA")
    f2_ax1.legend(labels=labels)
    

    f2_ax2 = fig2.add_subplot(gs[3:6, 0])
    # FA correlation dist
    f2_ax2 = sns.stripplot(x="TCK", y = "fa_y", data = correlation,
                        order = order, palette = new_palette,
                        rasterized=True, size=size)
    f2_ax2 = sns.boxplot(x="TCK", y = "fa_y", data = correlation,
                        order = order, palette = new_palette, showfliers=False
                        )
    
    f2_ax2.set(xlabel="fiber group label")
    f2_ax2.set(ylabel="r")
    f2_ax2.xaxis.set_major_locator(ticker.MultipleLocator(2))
    plt.xticks(rotation=45,ha='right')

    f2_ax3 = fig2.add_subplot(gs[0:2, 1])
    f2_ax3 = sns.stripplot(x="TCK", y = "bundle_adjacency_voxels", data = pairwise,
                        order= order, palette = new_palette,
                        rasterized=True, size=size)
    f2_ax3 = sns.boxplot(x="TCK", y = "bundle_adjacency_voxels", data = pairwise,
                        order= order, palette = new_palette,showfliers=False)   
    f2_ax3.set(xticklabels=[])
    f2_ax3.set(ylabel="bundle adjacency")
    f2_ax3.set(xlabel=None)

    #if btw=="compute": f2_ax3.set_yticks([0,0.5,1.0,1.5])

    f2_ax4 = fig2.add_subplot(gs[2:4, 1])
    f2_ax4 = sns.stripplot(x="TCK", y = "dice_voxels", data = pairwise,
                        order = order, palette = new_palette, 
                        rasterized=True, size=size)
    f2_ax4 = sns.boxplot(x="TCK", y = "dice_voxels", data = pairwise,
                        order = order, palette = new_palette, 
                        showfliers=False)
    f2_ax4.set(xticklabels=[])
    f2_ax4.set(ylabel="Dice index")
    f2_ax4.set(xlabel=None)

    f2_ax5 = fig2.add_subplot(gs[4:6, 1])
    f2_ax5 = sns.stripplot(x="TCK", y = "density_correlation", data = pairwise,
                        order = order,  palette = new_palette, 
                        rasterized=True, size=size)
    f2_ax5 = sns.boxplot(x="TCK", y = "density_correlation", data = pairwise,
                        order = order,  palette = new_palette, 
                        showfliers=False)
    f2_ax5.set(xlabel="fiber group label")   
    f2_ax5.set(ylabel="density correlation") 
    f2_ax5.xaxis.set_major_locator(ticker.MultipleLocator(2))
    plt.xticks(rotation=45,ha='right')

    return fig2 

groups = pd.read_csv(raw_csv / "groups.csv")
group = "IFG"
tck_to_plot = groups[group][groups[group].notna()]
profile_tmp = profile[profile["TCK"].isin(tck_to_plot)]
correlation_tmp = correlation[correlation["TCK"].isin(tck_to_plot)]                                                 
pairwise_tmp = pairwise[pairwise["TCK"].isin(tck_to_plot)]
palette = sns.color_palette("tab20")[0:10] + sns.color_palette("tab20b")

fig2 = plot_fig("compute", profile_tmp, correlation_tmp, pairwise_tmp, tck_to_plot)
#plt.show()
fig2.set_size_inches(21.5,13.24)

fig2.savefig( fig_dir / "Fig2_computational_IFG.svg", dpi=300, bbox_inches='tight')

fig2 = plot_fig("test-retest", profile_tmp, correlation_tmp, pairwise_tmp, tck_to_plot)
plt.show()
fig2.set_size_inches(9.88,5.93)
fig2.savefig( fig_dir / "Fig2_TRT_IFG.svg", dpi=300, bbox_inches='tight')


# generate shell script to visualize tract streamlines
# select only the left hemisphere

palette = palette[::2]
# palette = [matplotlib.colors.to_hex(x) for x in palette]

sub="S038"; ana="AL_07"
base_dir=f"/bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-{ana}/sub-{sub}/ses-T01/output"
tractparams = pd.read_csv(git_dir / "tractparams_AL_final_both_hemi.csv")
tract_dic = dict(zip(tractparams["slabel"], tractparams["roi2"]))


tcks = [x for i in tck_to_plot[::2] for x, y in tract_dic.items() if y==i]   
#tck_to_plot=tractparams.roi2
#tcks = [x for i in tck_to_plot[:47] for x, y in tract_dic.items() if y==i]   

# start writing shell script
with open("s5_visualization.sh", 'w') as f:
    f.write("#!/bin/bash\n")
    f.write(f"vglrun mrview \\\n")
    f.write(f"\t-load {base_dir}/flywheel/v0/output/RTP/fs/nans.nii.gz \\\n")
    f.write(f"\t-overlay.load {base_dir}/flywheel/v0/output/RTP/fs/brain_glass.nii.gz \\\n")
    f.write(f"\t-overlay.opacity 0.03 \\\n")
    f.write(f"\t-overlay.colourmap 0 \\\n")
    f.write(f"\t-overlay.load {base_dir}/flywheel/v0/output/RTP/fs/MD_outline.nii.gz \\\n")
    # f.write(f"\t-overlay.load {base_dir}/flywheel/v0/output/RTP/fs/ROIs/Left-MD_dil-1.nii.gz \\\n")    
    for tck, col, roi in zip(tcks, palette, tck_to_plot[::2]):
        col = tuple(x*255 for x in col)
        f.write(f"\t-overlay.load {base_dir}/flywheel/v0/output/RTP/fs/ROIs/{roi}_dil-1.nii.gz \\\n")
        f.write(f"\t-overlay.colour {col[0]},{col[1]},{col[2]} \\\n")
        f.write(f"\t-overlay.intensity 0.2,1 \\\n")
        f.write(f"\t-overlay.threshold_min 0.2 \\\n")
        f.write(f"\t-overlay.threshold_max 1 \\\n")
        f.write(f"\t-tractography.load {base_dir}/{tck}_clean.tck \\\n")
        f.write(f"\t-tractography.lighting 1 \\\n")
        f.write(f"\t-tractography.colour {col[0]},{col[1]},{col[2]} \\\n")
    f.write("\t-mode 3 -noannotations -fullscreen \n")
f.close()

# create colormap for surface ROI rendering
LUT_HCP = pd.read_csv(fig_dir / "LUT_HCP.txt", sep=" ")
LUT_HCP.columns = ["integral", "label"]
df_tmp = pd.DataFrame(columns = ["integral","label","R", "G", "B", "A"])
for group in groups.columns:
    if group=="AN": continue
    tck_to_plot = groups[group][groups[group].notna()]
    for col, roi in zip(palette, tck_to_plot[::2]):
        if roi=="L_Area_47l": roi="L_Area_47l_(47_lateral)"
        if len(LUT_HCP[LUT_HCP["label"]==roi]["integral"])==0:
            print(roi,col); continue
        integral = LUT_HCP[LUT_HCP["label"]==roi]["integral"].values[0]
        R = col[0]; G = col[1]; B = col[2]
        df_tmp = df_tmp.append({"integral":integral, "label":roi,
                         "R":R, "G":G, "B":B, "A":1},ignore_index = True)
df_tmp.to_csv(fig_dir / "HCP_surface_ROI_LUT.txt", sep="\t", index=False)
# cmd to attach this LUT to niml.dset
# MakeColorMap -usercolutfile HCP_surface_ROI_LUT.txt -suma_cmap toylut -sdset ~/suma_MNI152_2009/Orbital.niml.dset -overwrite


import tract_3D
import importlib
importlib.reload(tract_3D)
import tract_3D

x_cut = 128; y_cut = 118; z_cut = 153
camera = ((-705.8306541647723, 529.2247187019981, 425.76215289175894),
         (87.6887652403265, -14.697480338511443, 51.07560287191844), 
         (0.31958863044278235, -0.17523275927439252, 0.9312124287018455))
groups = pd.read_csv(raw_csv / "groups.csv")

with open('render_params.json') as d:
    render_params = json.load(d)

for group in groups.columns:
    tck_to_plot = groups[group][groups[group].notna()]
    tcks = {x:y for i in tck_to_plot[::2] for x, y in tract_dic.items() if y==i}
    for (key, value), col, n in zip(tcks.items(), palette, range(len(tcks))):
        view = tract_3D.tract_3D(bundle_filename=f"output/{key}_clean.tck", colors=col,
                            x_cut = x_cut, y_cut = y_cut, z_cut = z_cut,
                            camera = camera, interactive=True,
                            output=f"figures/{group}_{n:02d}_{value}.png")
        view["name"]=key; view["filename"]=value; view['group']=group
        view["color"]=col
        render_params[key]=view

import json
with open('render_params.json', 'w') as fp:
    json.dump(render_params, fp)

with open('render_params.json') as d:
    render_params = json.load(d)

for group in groups.columns:
    tck_to_plot = groups[group][groups[group].notna()]
    tcks = {x:y for i in tck_to_plot[::2] for x, y in tract_dic.items() if y==i}
    for (key, view), n in zip(render_params.items(), range(len(tcks))):
        x_cut = view["x_cut"]; y_cut = view["y_cut"]; z_cut = view["z_cut"]
        col = view["color"]; camera = view["camera"]
        tract_3D.tract_3D(bundle_filename=f"output/{key}_clean.tck", colors=col,
                            x_cut = x_cut, y_cut = y_cut, z_cut = z_cut,
                            camera = camera, interactive=False,
                            output=f"{group}_{n:02d}_{value}.png")


