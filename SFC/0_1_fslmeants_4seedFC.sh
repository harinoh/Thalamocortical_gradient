# !/bin/bash

ROI=PFCSS_Mot_cerebell68
opt=FWHM4_Parcel
p=1000
opt2=_BhatTHR

home_dir=/Volume/CCNC_T_4/harin_oh/1_thalamocortical
mask_dir=$home_dir/Gradient_result/Parcellation/3_1_Stat_anal/${ROI}_${opt}${p}${opt2}

for i in $@
do
	i_dir=$home_dir/SCZOHC/$i/REST/
	if [ ! -e $i_dir/${ROI}_L_thal_meants.txt ] ; then
		fslmeants -i $i_dir/rsfMRI_${i}_data.nii.gz -o $i_dir/${ROI}_L_thal_meants.txt \
			-m $mask_dir/ttest_brainstat_qvalue_G1_${ROI}_L_thal_ageCo.nii.gz
	else
		echo "subject ${i} meant time-series already calculated"
		continue
	fi
done
