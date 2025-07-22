#!/bin/bash

## If re-sampling to 3mm
group=SCZOHC
dir=/Volume/CCNC/harin_oh/1_thalamocortical/$group

#for i in N*
for i in $@
do
	rest_dir=$dir/$i/REST
	if [ ! -e $rest_dir/rsfMRI_${i}_prepFinal_FWHM4_3mm.nii.gz ] ; then
		flirt -in $rest_dir/rsfMRI_${i}_prepFinal_FWHM4.nii.gz -ref $rest_dir/rsfMRI_${i}_prepFinal_FWHM4.nii.gz -applyisoxfm 3 \
			-interp nearestneighbour -out $rest_dir/rsfMRI_${i}_prepFinal_FWHM4_3mm.nii.gz
	else
		echo "Subject ${i} already realigned"
		continue
	fi
done
