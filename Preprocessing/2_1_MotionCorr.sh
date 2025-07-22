# !/bin/bash

dir=/Volume/CCNC/harin_oh/1_thalamocortical

#for i in N*
for i in $@
do
        ## labelling directories ##
	
	REST=$dir/$i/REST
    Preproc=$dir/$i/preproc
    if [ ! -d $Preproc/MotOutlier ] ; then
        mkdir $Preproc/MotOutlier
    fi
	
	MOout=$Preproc/MotOutlier
	
	## 2. Motion Correction ##
	
	if [ ! -e $Preproc/rsfMRI_raw_MoCorr.nii.gz ] ; then
		mcflirt -in $Preproc/rsfMRI_raw.nii.gz -out $Preproc/rsfMRI_raw_MoCorr -mats -plots -rmsrel -rmsabs -spline_final
		echo "subject ${i} motion corrected"
	fi


	if [ ! -e $MOout/rsfMRI_raw_MOout_FDRMS_val ] ; then
		fsl_motion_outliers -i $Preproc/rsfMRI_raw.nii.gz -o $MOout/confound_MO -s $MOout/${i}_FDRMS_val_MO -p $MOout/${i}_FDRMS_1 --fdrms
		fsl_motion_outliers -i $Preproc/rsfMRI_raw_MoCorr.nii.gz -o $MOout/confound_MO_MoCorr -s $MOout/${i}_MoCorr_FDRMS_val_MO \
			-p $MOout/${i}_MoCorr_FDRMS_1 --fdrms
		echo "subject ${i} FDRMS calculated"
	fi

	echo "subject ${i} completed Motion correction (Step 2 of 9)"
done

