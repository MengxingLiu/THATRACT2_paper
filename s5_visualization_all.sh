#!/bin/bash
vglrun mrview \
	-load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/brainmask.nii.gz \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/Left-MD_dil-1.nii.gz \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_RetroSplenial_Complex_dil-1.nii.gz \
	-overlay.colour 31.0,119.0,180.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL1_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 31.0,119.0,180.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Entorhinal_Cortex_dil-1.nii.gz \
	-overlay.colour 255.0,127.0,14.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL2_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 255.0,127.0,14.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/Left-MammillaryBody_dil-1.nii.gz \
	-overlay.colour 44.0,160.0,44.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL3_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 44.0,160.0,44.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_TG_dorsal_dil-1.nii.gz \
	-overlay.colour 214.0,39.0,40.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL4_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 214.0,39.0,40.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_9_Middle_dil-1.nii.gz \
	-overlay.colour 148.0,103.0,189.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL5_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 148.0,103.0,189.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_10v_dil-1.nii.gz \
	-overlay.colour 140.0,86.0,75.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL6_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 140.0,86.0,75.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_p32_dil-1.nii.gz \
	-overlay.colour 227.0,119.0,194.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL7_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 227.0,119.0,194.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_10d_dil-1.nii.gz \
	-overlay.colour 127.0,127.0,127.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL8_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 127.0,127.0,127.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_anterior_10p_dil-1.nii.gz \
	-overlay.colour 188.0,189.0,34.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL9_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 188.0,189.0,34.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_11l_dil-1.nii.gz \
	-overlay.colour 23.0,190.0,207.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL10_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 23.0,190.0,207.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_47m_dil-1.nii.gz \
	-overlay.colour 57.0,59.0,121.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL11_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 57.0,59.0,121.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_13l_dil-1.nii.gz \
	-overlay.colour 107.0,110.0,207.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL12_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 107.0,110.0,207.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_posterior_OFC_Complex_dil-1.nii.gz \
	-overlay.colour 99.0,121.0,57.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL13_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 99.0,121.0,57.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_45_dil-1.nii.gz \
	-overlay.colour 181.0,207.0,107.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL14_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 181.0,207.0,107.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_posterior_47r_dil-1.nii.gz \
	-overlay.colour 140.0,109.0,49.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL15_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 140.0,109.0,49.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_IFSp_dil-1.nii.gz \
	-overlay.colour 231.0,186.0,82.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL16_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 231.0,186.0,82.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_IFJp_dil-1.nii.gz \
	-overlay.colour 132.0,60.0,57.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL17_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 132.0,60.0,57.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_9_Posterior_dil-1.nii.gz \
	-overlay.colour 214.0,97.0,107.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL18_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 214.0,97.0,107.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_9-46d_dil-1.nii.gz \
	-overlay.colour 123.0,65.0,115.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL19_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 123.0,65.0,115.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_46_dil-1.nii.gz \
	-overlay.colour 206.0,109.0,189.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL20_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 206.0,109.0,189.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Area_8C_dil-1.nii.gz \
	-overlay.colour 57.0,59.0,121.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL21_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 57.0,59.0,121.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Inferior_6-8_Transitional_Area_dil-1.nii.gz \
	-overlay.colour 107.0,110.0,207.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL22_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 107.0,110.0,207.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/L_Superior_Frontal_Language_Area_dil-1.nii.gz \
	-overlay.colour 99.0,121.0,57.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL23_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 99.0,121.0,57.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/R_RetroSplenial_Complex_dil-1.nii.gz \
	-overlay.colour 181.0,207.0,107.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL24_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 181.0,207.0,107.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/R_Entorhinal_Cortex_dil-1.nii.gz \
	-overlay.colour 140.0,109.0,49.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL25_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 140.0,109.0,49.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/Right-MammillaryBody_dil-1.nii.gz \
	-overlay.colour 231.0,186.0,82.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL26_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 231.0,186.0,82.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/R_Area_TG_dorsal_dil-1.nii.gz \
	-overlay.colour 132.0,60.0,57.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL27_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 132.0,60.0,57.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/R_Area_9_Middle_dil-1.nii.gz \
	-overlay.colour 214.0,97.0,107.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL28_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 214.0,97.0,107.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/R_Area_10v_dil-1.nii.gz \
	-overlay.colour 123.0,65.0,115.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL29_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 123.0,65.0,115.0 \
	-overlay.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/flywheel/v0/output/RTP/fs/ROIs/R_Area_p32_dil-1.nii.gz \
	-overlay.colour 206.0,109.0,189.0 \
	-overlay.intensity 0.2,1 \
	-overlay.threshold_min 0.2 \
	-overlay.threshold_max 1 \
	-tractography.load /bcbl/home/home_g-m/lmengxing/TESTDATA/analysis-AL_07/sub-S038/ses-T01/output/AL30_clean.tck \
	-tractography.lighting 1 \
	-tractography.colour 206.0,109.0,189.0 \
	-mode 3 -noannotations -fullscreen 