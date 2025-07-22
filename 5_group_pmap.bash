#!/bin/bash

group=SCZOHC
atlas=WholeCTX
sub_group=SCZ
opt=_FWHM4_IndivMul

## Defining directories
MAINDIR=/Volume/CCNC/harin_oh/1_thalamocortical
roitype='Gordon_L_thal_222 Gordon_R_thal_222'
OUTdir=$MAINDIR/Gradient_result/4_pmap
maskdir=$MAINDIR/ROI_atlas/mni_icbm152_nlin_sym_09b_nifti/ICBM152_nlin_sym_09b_gm_mask_${atlas}_flirt.nii.gz
pmap_dir=$MAINDIR/Gradient_result/4_pmap/${group}/${atlas}_${group}${opt}
output=$OUTdir/group/$sub_group/${atlas}_${sub_group}${opt}

if [ ! -d $output ] ; then
	mkdir -p $output
fi

mapfile -t subjects < $MAINDIR/SCZOHC/Subject_list_${sub_group}.txt

## 5-2. Group-level pmap construction
for j in $roitype
do
    if [ ! -e $OUTdir/Gordon_${j}_thal_222.pmaps.nii.gz ] ; then
        infiles=()
        for sub in "${subjects[@]}"
        do
            infiles+=("$DATADIR/${sub}_Gordon_${j}_thal_222_aligned_BhatTHR_mul.pmaps.nii.gz")
        done

        python $MAINDIR/code/congrads/pmap_from_Py_Aligned_cmaps_group.py \
            -i "${infiles[@]}" -d ${MAINDIR}/SCZOHC -g $sub_group \
            -o $OUTDIR/Gordon_${j}_thal_222
    else
        echo "ROI ${ROI} Group pmap alreayd created"
        continue
    fi
done
