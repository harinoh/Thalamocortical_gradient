#!/bin/bash

group=SCZOHC
atlas=Cerebellum
sub_group=SCZOHC
opt=_FWHM4

## Defining directories
dir=/Volume/CCNC_T_4/harin_oh/1_thalamocortical
roitype='Gordon_L_thal_222 Gordon_R_thal_222'
OUTdir=/Volume/CCNC_T_4/harin_oh/1_thalamocortical/Gradient_result/1_gradient_extraction
maskdir=$dir/ROI_atlas/mni_icbm152_nlin_sym_09b_nifti/ICBM152_nlin_sym_09b_gm_mask_${atlas}_flirt${opt2}.nii.gz

if [ ! -d $OUTdir/group/$sub_group/${atlas}${opt}${opt2} ] ; then
	mkdir -p $OUTdir/group/$sub_group/${atlas}${opt}${opt2}
fi


## 2. Group level gradient construction
for j in $roitype
do
	if [ ! -e $OUTdir/group/$sub_group/${atlas}${opt}/Gordon_R_thal_222.cmaps.nii.gz ] ; then
		python $dir/code/congrads/conmap_group_example.py \
		-r $dir/ROI_atlas/Gordon_Parcels/${j}.nii.gz \
       		-m $maskdir -o $OUTdir/group/$sub_group/${atlas}${opt}${opt2} --nmaps 10 --project \
		--group $group --sub_group $sub_group --opt ${opt}${opt2}
	else
		echo "${sub_group} group cmap already created"
		continue
	fi
done
