# !/bin/bash

dir=/Volume/CCNC_T_4/harin_oh/1_thalamocortical

#for i in N*
for i in $@
do
	## labelling directories ##

    Preproc=$dir/$i/preproc
    EPI=$Preproc/Bandpass/rsfMRI_${i}_bp.nii.gz

	Reg2MNI=$Preproc/Reg2MNI
	MNI_brain=/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz
	MNI_brain_mask=/usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz
#	ref=$dir/ROI_atlas/mni_icbm152_nlin_sym_09b_nifti/mni_icbm152_nlin_sym_09b/ICBM152_nlin_sym_09b_brain.nii.gz

	## 8-3. Registration to Standard space (MNI) _ epi_to_MNI ##

		# registering EPI image to MNI standard space
	if [ ! -e $Reg2MNI/epi2mni.nii.gz ] ; then
#		convertwarp --ref=$MNI_brain --warp1=$Reg2MNI/t1w2mni_warp --out=$Reg2MNI/epi2mni_warp #--premat=$Reg2MNI/t1w2mni_flirt.mat
		applywarp -i $EPI -r $MNI_brain -o $Reg2MNI/epi2mni -w $Reg2MNI/t1w2mni_warp
		
		echo "EPI to MNI completed"
	fi

	echo "subject ${i} moved to standard space"
done
