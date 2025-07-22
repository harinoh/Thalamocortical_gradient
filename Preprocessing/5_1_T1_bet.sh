#!/bin/bash

dir=/Volume/CCNC_T_4/harin_oh/1_thalamocortical

#for i in N*
for i in $@
do
	## labelling directories ##
	T1=$dir/$i/T1

	## 5-1. Co-registering functional image to structural image ##

	if [ ! -e $T1/T1_brain.nii.gz ] ; then
		bet $T1/T1.nii.gz $T1/T1_brain.nii.gz -f 0.3 -c 111 160 190 -m
	
	else
		echo "subject ${i} T1 brain extracted"
	fi

	if [ ! -e $T1/T1_brain_2mm_mask.nii.gz ] ; then
		flirt -in $T1/T1.nii.gz -out $T1/T1_2mm.nii.gz -ref $T1/T1.nii.gz -applyisoxfm 2
		flirt -in $T1/T1_brain.nii.gz -out $T1/T1_brain_2mm.nii.gz -ref $T1/T1_brain.nii.gz -applyisoxfm 2
		fslmaths $T1/T1_brain_2mm.nii.gz -bin $T1/T1_brain_2mm_mask.nii.gz
	fi

done

