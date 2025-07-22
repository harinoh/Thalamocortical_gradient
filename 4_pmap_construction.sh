#!/bin/bash

ROI=WholeCTX
opt=FWHM4_Parcel
p=1000
ref_group=SCZOHC

## Defining directories
MAINDIR=/Volume/CCNC_T_4/harin_oh/1_thalamocortical
DATADIR=$MAINDIR/SCZOHC
CMAP=$OUTDIR/Gradient_result/1_gradient_extraction/group/${ref_group}/${ROI}_${opt}${p}
L_CMAP=$CMAP/Gordon_L_thal_222_mul.cmaps.nii.gz
R_CMAP=$CMAP/Gordon_R_thal_222_mul.cmaps.nii.gz
OUTDIR=$OUTDIR/Gradient_result/4_pmap/SCZOHC/${ROI}_${opt}${p}

if [ ! -d $OUTDIR ] ; then
	mkdir -p $OUTDIR
fi


## 4. Constructing individual pmap using aligned cmap
for i in $@
do
	if [ ! -e $OUTDIR/${i}_Gordon_R_thal_222_aligned.pmaps.nii.gz ] ; then
		python $MAINDIR/code/congrads/pmap_from_Py_Aligned_cmaps_nifti.py \
			-i $DATADIR/$i/REST/rsfMRI_${i}_prepFinal_FWHM4.nii.gz \
			--cmap $L_CMAP \
			-m $MAINDIR/ROI_atlas/Schaefer2018_${p}Parcels_7Networks_order_FSLMNI152_2mm_${ROI}.nii.gz \
			-r $MAINDIR/ROI_atlas/Gordon_Parcels/Gordon_L_thal_222.nii.gz -o $OUTDIR --nmaps 5 --project

		python $MAINDIR/code/congrads/pmap_from_Py_Aligned_cmaps_nifti.py \
			-i $DATADIR/$i/REST/rsfMRI_${i}_prepFinal_FWHM4.nii.gz \
			--cmap $R_CMAP \
			-m $MAINDIR/ROI_atlas/Schaefer2018_${p}Parcels_7Networks_order_FSLMNI152_2mm_${ROI}.nii.gz \
			-r $MAINDIR/ROI_atlas/Gordon_Parcels/Gordon_R_thal_222.nii.gz -o $OUTDIR --nmaps 5 --project
	else
		echo "Subject ${i}'s pmaps already created"
		continue
	fi
done

