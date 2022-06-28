#!/bin/bash
vglrun mrview \
	-load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/brainmask.nii.gz \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL4_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 31.0,119.0,180.0 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL2_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 255.0,127.0,14.0 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL1_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 44.0,160.0,44.0 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL3_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 214.0,39.0,40.0 \
	-mode 3 -noannotations -fullscreen 
