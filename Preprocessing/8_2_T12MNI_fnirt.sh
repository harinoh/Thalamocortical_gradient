# !/bin/bash

dir=/Volume/CCNC_T_4/harin_oh/1_thalamocortical


#for i in N*
for i in $@
do
	## labelling directories ##

    Preproc=$dir/$i/preproc
	T1=$dir/$i/T1/T1_2mm.nii.gz
	T1_mask=$dir/$i/T1/T1_brain_2mm_mask.nii.gz
	T1_brain=$dir/$i/T1/T1_brain_2mm.nii.gz

	Reg2MNI=$Preproc/Reg2MNI
	MNI=/usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz
	MNI_brain=/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz
	MNI_brain_mask=/usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz
#	ref=$dir/ROI_atlas/mni_icbm152_nlin_sym_09b_nifti/mni_icbm152_nlin_sym_09b/ICBM152_nlin_sym_09b_brain.nii.gz

	## 8-2. Registration to Standard space (MNI) _ fnirt ##

		## deformation matrix from whole T1 -> whole MNI brain, applying deformation field to T1 brain extracted image
	if [ ! -e $Reg2MNI/t1w2mni_warp.nii.gz ] ; then
		fnirt --in=$T1 --inmask=$T1_mask --aff=$Reg2MNI/t1w2mni_flirt.mat \
			--cout=$Reg2MNI/t1w2mni_warp --jout=$Reg2MNI/t1w2mni_fnirt --iout=$Reg2MNI/t1w2mni_jac \
			--ref=$MNI --refmask=$MNI_brain_mask -v

		applywarp -i $T1_brain -r $MNI_brain -o $Reg2MNI/t1w2mni_brain_fnirt -w $Reg2MNI/t1w2mni_warp

		echo "T1 to MNI FNIRT completed in subject ${i}"
	fi

	echo "subject ${i} fnirt normalized to MNI, step 8-2 completed."
done
