# !/bin/bash

# Record time
start=$(date +%s)

h=L
atlas=WholeCTX
ref_group=SCZOHC
p=1000
opt=_FWHM4_Parcel${p}

## Defining directories
MAINDIR=/Volume/CCNC/harin_oh/1_thalamocortical
pmap_dir=$MAINDIR/Gradient_result
out_dir=$pmap_dir/5_pmap_alignment/SCZOHC/${atlas}${opt}
ref_pmap=$pmap_dir/4_pmap/group/$ref_group/${atlas}${opt}/Gordon_${h}_thal_222.pmaps.nii.gz
roiDIR=$MAINDIR/ROI_atlas/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm_${atlas}
if [ ! -e $out_dir ] ; then
        mkdir -p $out_dir
fi


## 6-2. Aligning individual pmap gradient to group pmap
if [ ! -e $out_dir/SPR9_Gordon_${h}_thal_222_aligned.pmaps.nii.gz ] ; then
	python $MAINDIR/code/congrads/alignment_thalamus_pmap_parcel.py \
		-sub $MAINDIR/SCZOHC/Subject_list_SCZOHC.txt \
		-roi ${roiDIR}.nii.gz -ref $ref_pmap \
		-pmap $pmap_dir/4_pmap/SCZOHC/${atlas}${opt} \
		-o $out_dir -opt $opt2 -nmaps 5
else
	echo "${atlas} pmap alignment is already complete"
fi

end=$(date +%s)
# Calculate difference
diff=$(( start - end ))
# Convert the difference in time
hours=$(( diff / 3600 ))
minutes=$(( (diff % 3600) / 60 ))
seconds=$(( diff % 60 ))
echo "For Subject $@ - Time taken: $hours hours, $minutes minutes, and $seconds seconds."
