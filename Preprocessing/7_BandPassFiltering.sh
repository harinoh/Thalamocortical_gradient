# !/bin/bash


dir=/Volume/CCNC/harin_oh/1_thalamocortical

#for i in N*
for i in $@
do
	## labelling directories #

    Preproc=$dir/$i/preproc

    if [ ! -d $Preproc/Bandpass ] ; then
        mkdir $Preproc/Bandpass
    fi

    Nuisance=$Preproc/NuiRe
    BP=$Preproc/Bandpass
        

	## 7. Band-pass filter ##
		#mean value across time as a baseline
	if [ ! -e $BP/tempMean.nii.gz ] ; then
		fslmaths $Nuisance/nuisance_regressed -Tmean $BP/tempMean
	fi

	if [ ! -e $BP/rsfMRI_${i}_bp.nii.gz ] ; then
		# sigma = (1 / Hz_cutoff) / ( 2 * TR)
			# 1 / Hz_cutoff = Hz_inSec (e.g. 0.1 Hz = 10 sec, 0.01 Hz = 100 sec) 
		fslmaths $Nuisance/nuisance_regressed -bptf 1.42857142857 14.2857142857 -add $BP/tempMean $BP/rsfMRI_${i}_bp
#		fslmaths ${BP}/rsfMRI_${i}_bp -Tmean ${BP}/mean_func_bp
	fi

	echo "Step 7 out of 9 preproc completed in subject ${i}"
done
