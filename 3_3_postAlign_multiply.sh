#!/bin/bash

ROI=WholeCTX
ref_grp=SCZOHC
p=1000

## Defining directories
homedir=/Volume/CCNC/harin_oh/1_thalamocortical
in_dir=$homedir/Gradient_result/2_Alignment/SCZOHC/${ROI}_${ref_grp}_FWHM4_Parcel${p}_IndivMul
thal_mask=$homedir/ROI_atlas/Gordon_Parcels
ROI_mask=$homedir/ROI_atlas/mni_icbm152_nlin_sym_09b_nifti/ICBM152_nlin_sym_09b_gm_mask_${ROI}_flirt.nii.gz

## Post alignment multiplication (to remove any small values outside the ROI)
for i in $@
do
	if [ ! -e $in_dir/${i}_Gordon_R_thal_222_aligned_mul.cmaps.nii.gz ] ; then
		fslmaths $in_dir/${i}_Gordon_L_thal_222_aligned.cmaps.nii.gz -mul $thal_mask/Gordon_L_thal_222.nii.gz \
			$in_dir/${i}_Gordon_L_thal_222_aligned_mul.cmaps.nii.gz

		fslmaths $in_dir/${i}_Gordon_R_thal_222_aligned.cmaps.nii.gz -mul $thal_mask/Gordon_R_thal_222.nii.gz \
			$in_dir/${i}_Gordon_R_thal_222_aligned_mul.cmaps.nii.gz
	else
		echo "subject ${i} cmap & pmap has been multiplied"
		continue
	fi
done
