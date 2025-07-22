#!/bin/bash

group=SCZOHC
atlas=WholeCTX
opt=_FWHM4
opt2=_Parcel
n=68
h=R

homedir=/Volume/CCNC/harin_oh/1_thalamocortical
OUTdir=$homedir/Gradient_result/1_gradient_extraction/$group/${atlas}${opt}${opt2}${n}
if [ ! -d $OUTdir ] ; then
        mkdir $OUTdir
fi

cd $homedir/$group

for ID in $@
do
    ## Defining directories
    REST=$homedir/$group/${ID}/REST/rsfMRI_${ID}_prepFinal${opt}.nii.gz
    inROI=$homedir/ROI_atlas/Gordon_Parcels/Gordon_${h}_thal_222.nii.gz
    outROI=$homedir/ROI_atlas/${atlas}${opt2}_flirt.nii.gz
#       outROI=$homedir/ROI_atlas/mni_icbm152_nlin_sym_09b_nifti/ICBM152_nlin_sym_09b_gm_mask_${atlas}_flirt${opt2}.nii.gz


    ## 1. Constructing individual gradient
    if [ ! -e $OUTdir/${ID}/Gordon_${h}_thal_222.cmaps.nii.gz ] ; then
        mkdir -p $OUTdir/$ID
        python $homedir/code/congrads/conmap_parcel_Bhat_alter.py -i ${REST} -r ${inROI} -m ${outROI} \
            -n $n -o $OUTdir/${ID} --nmaps 10 --project --roi_name $atlas

    else
        echo "subject ${ID} already done"
    fi

	if [ ! -e $OUTdir/${ID}/Gordon_${h}_thal_222_mul.cmaps.nii.gz ] ; then
		fslmaths $OUTdir/${ID}/Gordon_${h}_thal_222.cmaps.nii.gz -mul $inROI \
			$OUTdir/${ID}/Gordon_${h}_thal_222_mul.cmaps.nii.gz
		fslmaths $OUTdir/${ID}/Gordon_${h}_thal_222.pmaps.nii.gz -mul \
			$homedir/ROI_atlas/${atlas}${opt2}_flirt_mask.nii.gz \
			$OUTdir/${ID}/Gordon_${h}_thal_222_mul.pmaps.nii.gz
	else
		echo "ROI regions already multiplied"
		continue
	fi
done

