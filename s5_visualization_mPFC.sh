#!/bin/bash
vglrun mrview \
	-load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/nans.nii.gz \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/brain_glass.nii.gz \
	-overlay.opacity 0.03 \
	-overlay.colourmap 0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/Left-MDl_AND_Left-MDm.nii.gz \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_p32_dil-1.nii.gz \
	-overlay.colour 31.0,119.0,180.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL13_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 31.0,119.0,180.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_s32_dil-1.nii.gz \
	-overlay.colour 255.0,127.0,14.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL14_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 255.0,127.0,14.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_8BM_dil-1.nii.gz \
	-overlay.colour 44.0,160.0,44.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL8_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 44.0,160.0,44.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_9_Middle_dil-1.nii.gz \
	-overlay.colour 214.0,39.0,40.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL9_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 214.0,39.0,40.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_10v_dil-1.nii.gz \
	-overlay.colour 148.0,103.0,189.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL11_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 148.0,103.0,189.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_10r_dil-1.nii.gz \
	-overlay.colour 57.0,59.0,121.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL10_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 57.0,59.0,121.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_25_dil-1.nii.gz \
	-overlay.colour 107.0,110.0,207.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL12_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 107.0,110.0,207.0 \
	-mode 3 -noannotations -fullscreen 
