# !/bin/bash

MATLAB="/usr/local/MATLAB/R2020a/bin/matlab"
dir=/Volume/CCNC_T_4/harin_oh/1_thalamocortical/code/SFC_obesity

$MATLAB -nodisplay -nosplash -nodesktop -r "run('$dir/sfc_4_2.m');exit;"
