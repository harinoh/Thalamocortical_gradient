#!/bin/bash

ref_group=ECTOHC
ROI=WholeCTX
p=1000
opt=FWHM4_Parcel${p}


## Defining directories
homedir=/Volume/CCNC/harin_oh/1_thalamocortical
in_dir=$homedir/Gradient_result/1_gradient_extraction/group/$ref_group/${ROI}_${opt}
cmap_mask=$homedir/ROI_atlas/Gordon_Parcels
pmap_mask=$homedir/ROI_atlas/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm_${ROI}_mask.nii.gz

## Removing any small values (created from spyder)
if [ ! -e $in_dir/Gordon_Bi_thal_222_mul.cmaps.nii.gz ] ; then
	fslmaths $in_dir/Gordon_L_thal_222.cmaps.nii.gz -mul $cmap_mask/Gordon_L_thal_222.nii.gz \
		$in_dir/Gordon_L_thal_222_mul.cmaps.nii.gz
	fslmaths $in_dir/Gordon_R_thal_222.cmaps.nii.gz -mul $cmap_mask/Gordon_R_thal_222.nii.gz \
		$in_dir/Gordon_R_thal_222_mul.cmaps.nii.gz
	fslmaths $in_dir/Gordon_L_thal_222_mul.cmaps.nii.gz -add \
		$in_dir/Gordon_R_thal_222_mul.cmaps.nii.gz \
		$in_dir/Gordon_Bi_thal_222_mul.cmaps.nii.gz

        fslmaths $in_dir/Gordon_L_thal_222.pmaps.nii.gz -mul $pmap_mask \
                $in_dir/Gordon_L_thal_222_mul.pmaps.nii.gz
        fslmaths $in_dir/Gordon_R_thal_222.pmaps.nii.gz -mul $pmap_mask \
                $in_dir/Gordon_R_thal_222_mul.pmaps.nii.gz
else
        echo "group ${ref_group} using ROI ${ROI} pmap has been already sampled to surface"
fi
