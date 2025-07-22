%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Define seed regions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function seed_def()
basepath = '/Volume/CCNC/harin_oh/1_thalamocortical/';
maskpath = '/Volume/CCNC/harin_oh/1_thalamocortical/ROI_atlas/';
outpath = [basepath,'SFC_result/'];
if ~exist(outpath,'dir') % creating folder
    mkdir(outpath)
end

addpath(genpath('/Volume/CCNC/harin_oh/1_thalamocortical/code/SFC'));

%% 1) Define seed regions

% seed = MRIread([basepath,'Gradient_result/3_2_FunctionalNetwork/Thalamic_FN_Gordon_Bi_thal_222_SCZOHC_FWHM4.nii.gz']).vol;
% thal_seed = reshape(seed,902629,1);
% thal_seed = thal_seed(thal_comb_indx);
% 
% seed_idx_SMN = find(thal_seed==2);
% seed_idx_int = find(thal_seed==3 | thal_seed==4);
% seed_idx_HO = find(thal_seed==6 | thal_seed==7);
% 
% % save as mat file
% save([outpath,'FO_seed_region.mat'],'seed_idx_SMN');
% save([outpath,'Int_seed_region.mat'],'seed_idx_int');
% save([outpath,'HO_seed_region.mat'],'seed_idx_HO');

for Ctx_thal_comb_ver=1;
comb_seed = MRIread([maskpath,'mni_icbm152_nlin_sym_09b_nifti/ICBM152_nlin_sym_09b_gm_mask_WholeCTX_THAL_numlab_flirt.nii.gz']).vol;
thal_comb_seed = reshape(comb_seed,902629,1);

comb_seed_idx_whlThal = find(thal_comb_seed > 1);
comb_seed_idx_VIS = find(thal_comb_seed==2);
comb_seed_idx_SMN = find(thal_comb_seed==3);
comb_seed_idx_DAN = find(thal_comb_seed==4);
comb_seed_idx_VAN = find(thal_comb_seed==5);
comb_seed_idx_FPN = find(thal_comb_seed==7);
comb_seed_idx_DMN = find(thal_comb_seed==8);

comb_seed_idx_int = find(thal_comb_seed==4 | thal_comb_seed==5);
comb_seed_idx_HO = find(thal_comb_seed==7 | thal_comb_seed==8);

% save as mat file
save([outpath,'ctxthal_whlThal_seed_region.mat'],'comb_seed_idx_whlThal');
save([outpath,'ctxthal_VIS_seed_region.mat'],'comb_seed_idx_VIS');
save([outpath,'ctxthal_SMN_seed_region.mat'],'comb_seed_idx_SMN');
save([outpath,'ctxthal_DAN_seed_region.mat'],'comb_seed_idx_DAN');
save([outpath,'ctxthal_VAN_seed_region.mat'],'comb_seed_idx_VAN');
save([outpath,'ctxthal_FPN_seed_region.mat'],'comb_seed_idx_FPN');
save([outpath,'ctxthal_DMN_seed_region.mat'],'comb_seed_idx_DMN');

save([outpath,'ctxthal_Int_seed_region.mat'],'comb_seed_idx_int');
save([outpath,'ctxthal_HO_seed_region.mat'],'comb_seed_idx_HO');
end

% extracting seed index from thal_ctxthal
thal_mask = MRIread([maskpath,'Gordon_Parcels/Gordon_Bi_thal_222.nii.gz']).vol;
thal_mask = thal_mask > 0;
% ctx_mask = MRIread([maskpath, 'mni_icbm152_nlin_sym_09b_nifti/ICBM152_nlin_sym_09b_gm_mask_',num2str(roi),'_flirt.nii.gz']).vol;
% ctx_mask = ctx_mask > 0;
ctx_thal_mask = MRIread([maskpath,'mni_icbm152_nlin_sym_09b_nifti/ICBM152_nlin_sym_09b_gm_mask_WholeCTX_THAL_flirt.nii.gz']).vol;
ctx_thal_mask = ctx_thal_mask > 0;
indx_ctx_thal_mask = find(ctx_thal_mask > 0);

for x=1:size(thal_idx,1);
    idx = comb_seed(x);
    ctxthal_idx = find(indx_ctx_thal_mask == idx);
    indx_thal_ctxthal(x,1) = ctxthal_idx;
end
save([outpath,'thal_index_from_ctxthal.mat'],'indx_thal_ctxthal');

end
