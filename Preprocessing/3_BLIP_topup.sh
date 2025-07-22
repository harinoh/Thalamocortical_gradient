# !/bin/bash

dir=/Volume/CCNC/harin_oh/1_thalamocortical

for i in $@
do
	## Defining directories and files

	Preproc=$dir/$i/preproc
	BLIP_LR=$dir/$i/REST_LR/*nii.gz
	BLIP_RL=$dir/$i/REST_RL/*nii.gz
	REST=$Preproc/rsfMRI_raw_MoCorr

	if [ ! -d $Preproc/blip ] ; then
		mkdir $Preproc/blip
	fi
	
	BLIP=$Preproc/blip

	### 3. BLIP: Distortion correction

	if [ ! -e $BLIP/b0_blip_LR.nii.gz ] ; then
		fslroi $BLIP_LR $BLIP/b0_blip_LR 0 1 
		fslroi $BLIP_RL $BLIP/b0_blip_RL 0 1
		fslmerge -t $BLIP/both_b0 $BLIP/b0_blip_LR $BLIP/b0_blip_RL
	else
		echo "subject ${i} estimation field already created"
		
	fi

	if [ ! -e $BLIP/blip_correct_fieldcoef.nii.gz ] ; then
		topup --imain=$BLIP/both_b0 --datain=$dir/acq_param.txt --config=$dir/b02b0.cnf --out=$BLIP/blip_correct --iout=$BLIP/blip_unwarped.nii.gz
	else
		echo "subject ${i} topup finished"
	fi

	if [ ! -e $BLIP/rsfMRI_raw_Mocorr_blip.nii.gz ] ; then
		applytopup --imain=$REST --inindex=1 --datain=$dir/acq_param.txt --topup=$BLIP/blip_correct --method=jac --out=$BLIP/rsfMRI_raw_Mocorr_blip
		
	fi

	echo "subject ${i} completed in step 3 out of 9"
done
