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

