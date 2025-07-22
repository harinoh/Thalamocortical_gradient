#!/bin/bash

## Defining directories
homedir=/Volume/CCNC/harin_oh/1_thalamocortical
ROI=WholeCTX
in_dir=$homedir/Gradient_result/1_gradient_extraction/SCZOHC/${ROI}_FWHM4
thal_mask=$homedir/ROI_atlas/Gordon_Parcels
ROI_mask=$homedir/ROI_atlas/mni_icbm152_nlin_sym_09b_nifti/ICBM152_nlin_sym_09b_gm_mask_${ROI}_flirt.nii.gz

## Individual gradient also undergoing multiplication to remove any small values outside of ROI

for i in $@
do
	if [ ! -e $in_dir/$i/Gordon_R_thal_222_mul.cmaps.nii.gz ] ; then
		fslmaths $in_dir/$i/Gordon_L_thal_222.cmaps.nii.gz -mul $thal_mask/Gordon_L_thal_222.nii.gz \
			$in_dir/$i/Gordon_L_thal_222_mul.cmaps.nii.gz
		fslmaths $in_dir/$i/Gordon_R_thal_222.cmaps.nii.gz -mul $thal_mask/Gordon_R_thal_222.nii.gz \
			$in_dir/$i/Gordon_R_thal_222_mul.cmaps.nii.gz

		fslmaths $in_dir/$i/Gordon_L_thal_222.pmaps.nii.gz -mul $ROI_mask \
			$in_dir/$i/Gordon_L_thal_222_mul.pmaps.nii.gz
		fslmaths $in_dir/$i/Gordon_R_thal_222.pmaps.nii.gz -mul $ROI_mask \
			$in_dir/$i/Gordon_R_thal_222_mul.pmaps.nii.gz
	else
		echo "subject ${i} cmap & pmap has been multiplied"
		continue
	fi
done
