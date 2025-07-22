# !/bin/bash

dir=/Volume/CCNC/harin_oh/1_thalamocortical

#for i in E*
for i in $@
do
        ## labelling directories ##

	Preproc=$dir/$i/preproc
	REST=$Preproc/blip/rsfMRI_raw_Mocorr_blip
	OUT=$Preproc/rsfMRI_raw_MoCorr_blip_bet

	## 4. Brain extraction ##

	if [ ! -e ${OUT}.nii.gz ] ; then
		fslmaths $REST -Tmean $Preproc/mean_func
		bet2 $Preproc/mean_func $Preproc/mask -f 0.4 -n -m 
		fslmaths $REST -mas $Preproc/mask_mask $OUT
	else
		echo "subject ${i} fMRI bet already completed"
		continue
	fi
	
	echo "subject ${i} completed step 4 of 9"

done
