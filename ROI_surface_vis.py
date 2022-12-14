
import subprocess as sp
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


groups = pd.read_csv(raw_csv / "groups.csv")

suma_dir = '/home/mengxing/suma_MNI152_2009'
LUT = f"{suma_dir}/LUT_HCP.txt"

HCP_surf = f"{suma_dir}/lh.std.141.Glasser_HCP.lbl.niml.dset"

LUT = pd.read_csv(LUT, delim_whitespace=True)
LUT.columns = ["value", "label"]

for col in groups.columns:

    ind = []
    for i in groups[col]:
        #print(i)
        if pd.isna(i):
            continue
        if i[0:2] == "L_":
           
            value = LUT[LUT["label"]==i]["value"].values[0]
            ind.append(value)

    ind = ','.join(str(n) for n in ind)
    cmd_str = f"3dcalc -a {HCP_surf} -expr 'a * amongst(a,{ind})' -prefix {suma_dir}/{col}.niml.dset -overwrite"
    print(cmd_str)
    sp.call(cmd_str,shell=True)
    cmd_str = f"3drefit -copyaux {HCP_surf} {suma_dir}/{col}.niml.dset -overwrite"
    print(cmd_str)
    sp.call(cmd_str,shell=True)


