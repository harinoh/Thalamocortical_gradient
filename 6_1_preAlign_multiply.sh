#!/bin/bash

ref_group=SCZOHC
ROI=WholeCTX
p=1000
opt=FWHM4_Parcel${p}

## Defining directories
homedir=/Volume/CCNC/harin_oh/1_thalamocortical
in_dir=$homedir/Gradient_result/4_pmap/SCZOHC/${ROI}_${opt}
ROI_mask=$homedir/ROI_atlas/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm_${ROI}_mask.nii.gz

## 6-1. Removing small values outside of ROI
for i in $@
do
	if [ ! -e $in_dir/${i}_Gordon_R_thal_222_aligned_mul.pmaps.nii.gz ] ; then
		fslmaths $in_dir/${i}_Gordon_L_thal_222_aligned.pmaps.nii.gz -mul $ROI_mask \
			$in_dir/${i}_Gordon_L_thal_222_aligned_mul.pmaps.nii.gz
		fslmaths $in_dir/${i}_Gordon_R_thal_222_aligned.pmaps.nii.gz -mul $ROI_mask \
			$in_dir/${i}_Gordon_R_thal_222_aligned_mul.pmaps.nii.gz
	else
		echo "subject ${i} pmap has been multiplied"
		continue
	fi
done
