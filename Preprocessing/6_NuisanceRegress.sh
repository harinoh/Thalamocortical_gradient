#!/bin/bash

dir=/Volume/CCNC/harin_oh/1_thalamocortical


#for i in S*
for i in $@
do
	## labelling directories ##

    Preproc=$dir/$i/preproc

    if [ ! -d $Preproc/NuiRe ] ; then
        mkdir $Preproc/NuiRe
    fi
	
	Nuisance=$Preproc/NuiRe
	EPI_2T1=$Preproc/epi2t1w.nii.gz

	## 6. Nuisance Signal regression ##

	if [ ! -e $Nuisance/nuisance.mat ] ; then
		flirt -in $Preproc/epi2t1w_fast_pve_0.nii.gz -ref $EPI_2T1 -applyxfm -init $Preproc/t1w2epi.mat -out $Nuisance/csf2epi.nii.gz
		flirt -in $Preproc/epi2t1w_fast_pve_2.nii.gz -ref $EPI_2T1 -applyxfm -init $Preproc/t1w2epi.mat -out $Nuisance/wm2epi.nii.gz
		fslmaths $Nuisance/csf2epi.nii.gz -thr 0.5 $Nuisance/50thr_csf2epi.nii.gz
		fslmaths $Nuisance/wm2epi.nii.gz -thr 0.75 $Nuisance/75thr_wm2epi.nii.gz
		fslmeants -i $EPI_2T1 -o $Nuisance/CSF_noise.txt -m $Nuisance/50thr_csf2epi.nii.gz
		fslmeants -i $EPI_2T1 -o $Nuisance/WM_noise.txt -m $Nuisance/75thr_wm2epi.nii.gz

			##collating all motion correction parameters
		paste $Preproc/rsfMRI_raw_MoCorr.par $Nuisance/CSF_noise.txt $Nuisance/WM_noise.txt > $Nuisance/nuisance.txt
		Text2Vest $Nuisance/nuisance.txt $Nuisance/nuisance.mat
	else
		echo "nuisance.mat already exist in ${i}"
	fi

		## regressing out CSF, wm and motion correction
	if [ ! -e $Nuisance/nuisance_regressed.nii.gz ] ; then
		fsl_regfilt -i $EPI_2T1 -d $Nuisance/nuisance.mat -o $Nuisance/nuisance_regressed -f 1,2,3,4,5,6,7,8
#		fsl_regfilt -i ${EPI_Standard} -d ${Nuisance}/nuisance.mat -o $Nuisance/nuisance_regressed_derivativesOnly -f 1,2,3
	else
		echo "nuisance_regression completed in ${i}"
		continue
	fi

	echo "Step 6 of preproc completed in subject ${i}"
done
