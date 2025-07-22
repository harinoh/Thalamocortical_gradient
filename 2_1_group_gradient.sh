#!/bin/bash

atlas=WholeCTX
sub_group=SCZOHC
p=1000
opt=_FWHM4_Parcel${p}

## Defining directories
homedir=/Volume/CCNC/harin_oh/1_thalamocortical
roitype='Gordon_L_thal_222 Gordon_R_thal_222'
OUTdir=$homedir/Gradient_result/1_gradient_extraction/group/$sub_group/${atlas}${opt}
maskdir=$homedir/ROI_atlas/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm_${atlas}.nii.gz

if [ ! -d $OUTdir ] ; then
	mkdir -p $OUTdir
fi

mapfile -t subjects < $homedir/SCZOHC/Subject_list_${sub_group}.txt


## 2. Group-level gradient construction
for j in $roitype
do
	if [ ! -e $OUTdir/Gordon_R_thal_222_Bhat.cmaps.nii.gz ] ; then
		infiles=()
                for sub in "${subjects[@]}"
                do
                        infiles+=("$homedir/SCZOHC/${sub}/REST/rsfMRI_${sub}_prepFinal_FWHM4.nii.gz")
                done

		python $homedir/code/congrads/conmap_group_Bhat.py \
		-r $homedir/ROI_atlas/Gordon_Parcels/${j}.nii.gz \
       		-i "${infiles[@]}" \
		-m $maskdir -o $OUTdir --nmaps 10 --project
	else
		echo "${sub_group} group cmap already created"
		continue
	fi
done
