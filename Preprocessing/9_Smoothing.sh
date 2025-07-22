#!/bin/bash

dir=/Volume/CCNC_T_4/harin_oh/1_thalamocortical

for i in $@
do
	## labeling directories ##

	EPI2MNI=$dir/$i/preproc/Reg2MNI/epi2mni.nii.gz
	REST=$dir/$i/REST

	## 9. Smoothing by 6mm FWHM ##

	if [ ! -e $REST/rsfMRI_${i}_prepFinal_FWHM4.nii.gz ] ; then
		# sigma = wanted FWHM / 2.354
		
		#fslmaths $EPI2MNI -s 2.55 $REST/rsfMRI_${i}_prepFinal.nii.gz
		fslmaths $EPI2MNI -s 1.7 $REST/rsfMRI_${i}_prepFinal_FWHM4.nii.gz
	fi

	echo "Final step of preproc complated in subject ${i}"
done
