#  !/bin/bash

dir=/Volume/CCNC/harin_oh/1_thalamocortical

#for i in N*
for i in $@
do
	## labelling directories ##
	T1=$dir/$i/T1/*2*.nii.gz
	REST=$dir/$i/REST
	
	Preproc=$dir/$i/preproc
	if [ ! -d $Preproc ] ; then
		mkdir $Preproc
	fi

	if [ ! -e $dir/$i/REST/REST.nii.gz ] ; then
		cp $REST/*2*nii.gz $REST/REST.nii.gz
		cp $T1 $dir/$i/T1/T1.nii.gz
	fi

	raw=$REST/REST.nii.gz

	## 1.Discard initial 4 ts ##

	if [ ! -e $Preproc/rsfMRI_raw.nii.gz ] ; then
		fslroi $raw $Preproc/rsfMRI_raw.nii.gz 4 246
	fi

	echo "subject ${i} completed Step 1 of 9"
done
