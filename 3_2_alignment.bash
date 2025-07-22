#!/bin/bash

atlas=WholeCTX
ref_group=SCZOHC
opt=_FWHM4_Parcel
h=L
p=1000

## Defining directories
MAINDIR='/Volume/CCNC_T_4/harin_oh/1_thalamocortical'
INDIR=$MAINDIR/Gradient_result/1_gradient_extraction/SCZOHC/${atlas}${opt}${p}
OUTDIR=$MAINDIR/Gradient_result/2_Alignment/SCZOHC/${atlas}_${ref_group}${opt}${p}_IndivMul
refIMG=$MAINDIR/Gradient_result/1_gradient_extraction/group/$ref_group/${atlas}${opt}${p}
inROI=$MAINDIR/ROI_atlas/Gordon_Parcels/Gordon_${h}_thal_222.nii.gz

if [ ! -d $OUTDIR ] ; then
    mkdir -p $OUTDIR
fi


## 3. Aligning individual gradient to the group gradient
if [ ! -e $OUTDIR/SPR9_Gordon_${h}_thal_222_aligned.cmaps.nii.gz ] ; then
	python $MAINDIR/code/congrads/alignment_thalamus.py -sub $MAINDIR/SCZOHC/Subject_list_${ref_group}.txt \
		-roi $inROI -ref $refIMG/Gordon_${h}_thal_222_mul.cmaps.nii.gz \
		-o $OUTDIR --indivCMAP $INDIR
else
	echo "all subject for ${atlas} has already been aligned"
	continue
fi

