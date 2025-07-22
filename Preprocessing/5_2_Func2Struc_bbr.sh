#!/bin/bash


dir=/Volume/CCNC/harin_oh/1_thalamocortical


#for i in SPR*
for i in $@
do
	## labelling directories ##
	echo "performing epi_reg on ${i}"

    Preproc=$dir/$i/preproc

	if [ -e $Preproc/rsfMRI_raw_MoCorr_blip_bet.nii.gz ] ; then
		EPI=$Preproc/rsfMRI_raw_MoCorr_blip_bet

	else
		echo "subject ${i} did not perform blip"
		EPI=$Preproc/rsfMRI_raw_MoCorr_bet
	fi

	T1=$dir/$i/T1

	## 5-2. Co-registering functional image to structural image ##

	if [ ! -e $Preproc/epi2t1w.nii.gz ] ; then
		epi_reg --noclean --epi=$EPI --t1=$T1/T1_2mm --t1brain=$T1/T1_brain_2mm --out=$Preproc/epi2t1w	
	fi
	
	if [ ! -e ${Preproc}/t1w2epi.mat ] ; then
	## for future when transforming MNI standard space transformation to functional image
		convert_xfm -inverse -omat ${Preproc}/t1w2epi.mat ${Preproc}/epi2t1w.mat 
		
	else
		echo "check for inputs for epi_reg in subject ${i}"
		continue
	fi

	echo "Step 5 of Preproc, completed in subject ${i}"
done

