# !/bin/bash

MATLAB="/usr/local/MATLAB/R2020a/bin/matlab"
dir=/Volume/CCNC/harin_oh/1_thalamocortical/code/SFC

$MATLAB -nodisplay -nosplash -nodesktop -r "run('$dir/binconn_1.m');exit;"
