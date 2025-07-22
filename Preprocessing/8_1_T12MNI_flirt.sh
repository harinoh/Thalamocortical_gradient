# !/bin/bash

dir=/Volume/CCNC_T_4/harin_oh/1_thalamocortical

#for i in N*
for i in $@
do
	## labelling directories ##
    Preproc=$dir/$i/preproc
	T1_brain=$dir/$i/T1/T1_brain_2mm
	EPI2T1w=${Preproc}/epi2t1w.mat
	
	if [ ! -d $Preproc/Reg2MNI ] ; then
		mkdir $Preproc/Reg2MNI
	fi

	Reg2MNI=$Preproc/Reg2MNI
	MNI_brain=/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz
	MNI_brain_mask=/usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz


	## 8. Registration to Standard space (MNI) ##

	if [ ! -e $Reg2MNI/t1w2mni_flirt.nii.gz ] ; then
		flirt -in $T1_brain -ref $MNI_brain -out $Reg2MNI/t1w2mni_flirt -omat $Reg2MNI/t1w2mni_flirt.mat -searchrx -180 180 \
			#-searchrx -180 180 -searchrz -180 180 -v		# searchrx/y/z for orientation
		
		convert_xfm -inverse -omat $Reg2MNI/mni2t1w_flirt.mat $Reg2MNI/t1w2mni_flirt.mat
		
		echo "T1 to MNI FLIRT completed"
	fi

		##concatenating transformation matrix of EPI -> T1 & T1 -> MNI
	if [ ! -e $Reg2MNI/epi2mni.mat ] ; then
		convert_xfm -omat $Reg2MNI/epi2mni.mat -concat $Reg2MNI/t1w2mni_flirt.mat $EPI2T1w
		convert_xfm -inverse -omat $Reg2MNI/mni2epi.mat $Reg2MNI/epi2mni.mat
	fi

	echo "subject ${i} completed preproc step 8 out of 9"
done
