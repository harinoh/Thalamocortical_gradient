# !/bin/bash

dir=/Volume/CCNC/harin_oh/1_thalamocortical

cd $dir/QC/2_MotionCorrect

for i in *_motion_outliers.txt
do
	echo "${i}"
	fgrep -o 0 ${i} | wc -l
done
