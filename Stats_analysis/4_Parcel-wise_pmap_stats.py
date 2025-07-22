#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 08:52:08 2024

@author: harin_oh
"""


import nibabel as nib
import numpy as np
import pandas as pd
import os
import sys

ROI= 'WholeCTX'
ref_group='SCZOHC'
par_num=1000
opt=f'_FWHM4_Parcel{par_num}' 
sub_group1 = 'OHC' # NOR
sub_group2 = 'SCZ' # patient
home_dir = '/Volume/CCNC/harin_oh/1_thalamocortical'
pmap_path = f'{home_dir}/Gradient_result/5_pmap_alignment/SCZOHC/{ROI}{opt}'
file_name_group1 = f'{home_dir}/SCZOHC/Subject_list_{sub_group1}.txt'
file_name_group2 = f'{home_dir}/SCZOHC/Subject_list_{sub_group2}.txt'
out_dir = f'{home_dir}/Gradient_result/6_1_Stats/SCZOHC/{ROI}{opt}'

brainstat_path = os.path.abspath(os.path.join(f'{home_dir}/code/BrainStat-master'))

if brainstat_path not in sys.path:
    sys.path.append(brainstat_path)

def read_subject_list(file_name):
    with open(file_name, 'r') as my_file:
        file_handle = my_file.read()
        subjects = file_handle.rstrip()
        subject_list = subjects.split("\n")
    return subject_list

def nii_read(subject_list):
    G1_L_thal = {}
    G2_L_thal = {}
    G3_L_thal = {}
    
    G1_R_thal = {}
    G2_R_thal = {}
    G3_R_thal = {}
    for sub in subject_list:
        L_thal = nib.load(f'{pmap_path}/{sub}_Gordon_L_thal_222_aligned_mul.pmaps.nii.gz')
        R_thal = nib.load(f'{pmap_path}/{sub}_Gordon_R_thal_222_aligned_mul.pmaps.nii.gz')
        
        L_thal_G1 = L_thal.get_fdata()[:, :, :, 0]
        R_thal_G1 = R_thal.get_fdata()[:, :, :, 0]
        L_thal_G2 = L_thal.get_fdata()[:, :, :, 1]
        R_thal_G2 = R_thal.get_fdata()[:, :, :, 1]
        L_thal_G3 = L_thal.get_fdata()[:, :, :, 2]
        R_thal_G3 = R_thal.get_fdata()[:, :, :, 2]
        
        G1_L_thal[sub] = L_thal_G1
        G2_L_thal[sub] = L_thal_G2
        G3_L_thal[sub] = L_thal_G3
        
        G1_R_thal[sub] = R_thal_G1
        G2_R_thal[sub] = R_thal_G2
        G3_R_thal[sub] = R_thal_G3
        
    return G1_L_thal, G2_L_thal, G1_R_thal, G2_R_thal, G3_L_thal, G3_R_thal


#### Main ###

if not os.path.exists(f'{out_dir}/Figure'):
    os.makedirs(f'{out_dir}/Figure')
    print(f'Folder {out_dir}/Figure created.')
else:
    print(f'Folder {out_dir}/Figure already exists')
    
    
print('Loading fMRI images & thalamic atlas')
# read subject file
subject_list_group1 = read_subject_list(file_name_group1)
subject_list_group2 = read_subject_list(file_name_group2)
print(f'NOR group: {subject_list_group1}')
print(f'patient group: {subject_list_group2}')

# Loading thalamic gradient
G1_L_data1, G2_L_data1, G3_L_data1, G1_R_data1, G2_R_data1, G3_R_data1 = nii_read(subject_list_group1)      # group NOR
G1_L_data2, G2_L_data2, G3_L_data2, G1_R_data2, G2_R_data2, G3_R_data2 = nii_read(subject_list_group2)      # group FEP

# Loading demographic data
cols = ['Folder_name', 'Group','age', 'sex_num'] #, 'sex_num'
demo = pd.read_excel(f'{home_dir}/BCS_ECT_OHC_rsfMRI_Demo.xlsx', sheet_name='baseline_n129', usecols=cols, index_col='Folder_name')

# Load roi & thalamic atlas
roi_atlas = nib.load(f'{home_dir}/ROI_atlas/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm_{ROI}.nii.gz')
roi_data = roi_atlas.get_fdata()
roidims = roi_data.shape
nVoxels = np.prod(roidims)
roi_data = np.reshape(roi_data,nVoxels)

par_val, par_counts = np.unique(roi_data, return_counts=True) 
par_val = par_val[par_val != 0]

    
    
###############################
#### 1. Extract cortical region
def roi_extract(G_data, roi_data, subject_list):
    G_data_roi={}
    G_data_roi_list = []
    for sub in subject_list:
        sub_data = G_data[sub]
        sub_data = np.reshape(sub_data, nVoxels)
        mean_par = np.zeros(len(par_val))
        for a, p in enumerate(par_val):
            roi_indx = np.where(roi_data == p)[0]
            sub_roi = sub_data[roi_indx]
            mean_par[a] = np.mean(sub_roi)
        G_data_roi[sub] = mean_par
        G_data_roi_list.append(mean_par)
    return G_data_roi, G_data_roi_list

# group1
G1_L_data1_pmap, G1_L_data1_pmap_list = roi_extract(G1_L_data1, roi_data, subject_list_group1)
G1_R_data1_pmap, G1_R_data1_pmap_list = roi_extract(G1_R_data1, roi_data, subject_list_group1)

# group2
G1_L_data2_pmap, G1_L_data2_pmap_list = roi_extract(G1_L_data2, roi_data, subject_list_group2)
G1_R_data2_pmap, G1_R_data2_pmap_list = roi_extract(G1_R_data2, roi_data, subject_list_group2)

del G1_L_data1, G1_R_data1, G1_L_data2, G1_R_data2


#### 1-2. Normalization
def norm(G_data,subject_list):
    G_data_norm = {}
    G_data_norm_list = []
    for sub in subject_list:
        sub_data = G_data[sub]
        sub_min = np.min(sub_data)
        sub_max = np.max(sub_data)
        
        norm_data = 2 * (sub_data - sub_min) / (sub_max - sub_min) - 1
        G_data_norm_list.append(norm_data)
        G_data_norm[sub] = norm_data
    
    return G_data_norm, G_data_norm_list
# group1
G1_L_data1_pmap_norm, G1_L_data1_pmap_norm_list = norm(G1_L_data1_pmap, subject_list_group1)
G1_R_data1_pmap_norm, G1_R_data1_pmap_norm_list = norm(G1_R_data1_pmap, subject_list_group1)

# group2
G1_L_data2_pmap_norm, G1_L_data2_pmap_norm_list = norm(G1_L_data2_pmap, subject_list_group2)
G1_R_data2_pmap_norm, G1_R_data2_pmap_norm_list = norm(G1_R_data2_pmap, subject_list_group2)



###############################
#### 2. Regression
# Regress out
from sklearn.linear_model import LinearRegression

def regress_glm(G_L_data1, G_L_data2):
    G1_data1_arr = np.array(list(G_L_data1.values()))
    G1_data2_arr = np.array(list(G_L_data2.values()))
    all_array = np.row_stack((G1_data1_arr, G1_data2_arr))
    
    # Covariate + 1 designated to X
    X = np.column_stack((np.ones_like(demo['age']), demo['Group'], demo['age'])) #  , demo['sex_num']
    model = LinearRegression()
    model.fit(X, all_array)
    
    resid = all_array - model.predict(X)
    data1_arr = resid[:len(G_L_data1),:]
    data2_arr = resid[len(G_L_data1):,:]
    
    data1_dic = {k: data1_arr[i] for i, k in enumerate(G_L_data1)}
    data2_dic = {k: data2_arr[i] for i, k in enumerate(G_L_data2)}
    
    return resid, data1_dic, data2_dic

G1_resid, G1_L_data1_resid, G1_L_data2_resid = regress_glm(G1_L_data1_pmap, G1_L_data2_pmap)
G1_resid, G1_R_data1_resid, G1_R_data2_resid = regress_glm(G1_R_data1_pmap, G1_R_data2_pmap)



###############################
### 3. Statistical testing
from scipy import stats
from statsmodels.stats.multitest import multipletests

## t-test version
def ttest_construct(G_L_data1, G_L_data2):
    G_t_stat = []
    G_p_val = []
    for vox in range(0,len(par_val)):
        G_L_data1_arr = np.array([v[vox] for k, v in G_L_data1.items()])
        G_L_data2_arr = np.array([v[vox] for k, v in G_L_data2.items()])
        
        t_stat, pval = stats.ttest_ind(G_L_data1_arr, G_L_data2_arr)
        G_t_stat.append(t_stat)
        G_p_val.append(pval)
        
    uncorr_p = np.where(np.array(G_p_val) < 0.05)[0]
    print('uncorrected voxels:', len(uncorr_p))
    FDR_corr_p_val = multipletests(G_p_val, method='fdr_bh', alpha = 0.05)[1]
    FDR_corr_p = np.where(FDR_corr_p_val < 0.05)[0]
    print('FDR_corr voxels:', len(FDR_corr_p))
    
    return G_t_stat, G_p_val, uncorr_p, FDR_corr_p

# w/o regression
G1_L_t_stat, G1_L_p_val, uncorr_p_G1_L, FDR_corr_p_G1_L = ttest_construct(G1_L_data1_pmap, G1_L_data2_pmap)
G1_R_t_stat, G1_R_p_val, uncorr_p_G1_R, FDR_corr_p_G1_R = ttest_construct(G1_R_data1_pmap, G1_R_data2_pmap)

# w/ regression
G1_L_t_stat_ageCo, G1_L_p_val_ageCo, uncorr_p_G1_L_ageCo, FDR_corr_p_G1_L_ageCo = ttest_construct(G1_L_data1_resid, G1_L_data2_resid)
G1_R_t_stat_ageCo, G1_R_p_val_ageCo, uncorr_p_G1_R_ageCo, FDR_corr_p_G1_R_ageCo = ttest_construct(G1_R_data1_resid, G1_R_data2_resid)



## Surfstat version
from brainstat.stats.SLM import SLM
from brainstat.stats import terms

def brainstat_ttest(G_L_data1_thal_list, G_L_data2_thal_list):
    term_intercept = terms.FixedEffect(1, names='intercept')
    term_age = terms.FixedEffect(demo['age'], names='age')
    term_group = terms.FixedEffect(demo['Group'], names='group')
    term_sex = terms.FixedEffect(demo['sex_num'], names='sex')
    #term_fdrms = terms.FixedEffect(concat_mean_fdrms)
    
    group = demo['Group']
    model = term_intercept + term_group + term_age #+ term_sex # + term_fdrms 
    slm_age = SLM(model, -group, correction='fdr')

    G_L_concat = np.array(G_L_data1_thal_list + G_L_data2_thal_list)

    slm_age.fit(G_L_concat)
    ttest = slm_age.t
    qval = slm_age.Q
    fdr_pval_idx = np.where(qval < 0.05)[0]
    
    if slm_age.P is None:
        print("P-values are not calculated by default. Attempting manual calculation...")

        # Assuming T-statistics are available in the fitted model
        if hasattr(slm_age, 't'):
            from scipy.stats import t

            df = slm_age.df  # Degrees of freedom
            t_values = slm_age.t   # T-statistics

            # Calculate p-values
            p_values = 2 * (1 - t.cdf(np.abs(t_values), df))
            slm_age.P = p_values
            pval = np.transpose(p_values)
        else:
            print("T-statistics not available. Cannot calculate P-values manually.")
    else:
        print("P-values are calculated.")
    
    print('uncorrected p-value:',len(np.where(pval < 0.05)[0]))
    print('FDR corrected p-value:',len(np.where(qval < 0.05)[0]))
    
    return ttest, qval, fdr_pval_idx, pval

G1_L_ttest, G1_L_qval, G1_L_qval_idx, G1_L_pval = brainstat_ttest(G1_L_data1_pmap_list, G1_L_data2_pmap_list)
G1_R_ttest, G1_R_qval, G1_R_qval_idx, G1_R_pval = brainstat_ttest(G1_R_data1_pmap_list, G1_R_data2_pmap_list)

# normalization version
G1_L_ttest, G1_L_qval, G1_L_qval_idx, G1_L_pval = brainstat_ttest(G1_L_data1_pmap_norm_list, G1_L_data2_pmap_norm_list)
G1_R_ttest, G1_R_qval, G1_R_qval_idx, G1_R_pval = brainstat_ttest(G1_R_data1_pmap_norm_list, G1_R_data2_pmap_norm_list)



###############################
## 4. rebuild into 3D image to save
roi_affine = roi_atlas.affine

def reconstruct_mask_2nifti(qval_idx):
    mask = np.zeros(len(par_val))
    mask[qval_idx] = 1

    mask_3d = np.zeros(roi_atlas.shape)
    mask_3d = np.reshape(mask_3d,nVoxels)
    for a, p in enumerate(par_val):
        mask_par = mask[a]
        ind = np.where(roi_data == p)[0]
        mask_3d[ind] = mask_par
    mask_3d = np.reshape(mask_3d, roi_atlas.shape)
    
    return mask_3d


# saving ttest result
mask_G1_L_3d = reconstruct_mask_2nifti(G1_L_qval_idx)
mask_G1_R_3d = reconstruct_mask_2nifti(G1_R_qval_idx)

def save_2nifti(niftiImg,outname):
    nifti_img = nib.Nifti1Image(niftiImg, roi_affine)
    nib.save(nifti_img,os.path.join(f'{out_dir}/{outname}.nii.gz'))

# L-thalamus
print(out_dir)
opt2='_noRegress_norm'
save_2nifti(mask_G1_L_3d,f'ttest_G1_{ROI}_L_significantDiff_brainstat{opt2}')
save_2nifti(mask_G2_L_3d,f'ttest_G2_{ROI}_L_significantDiff_brainstat{opt2}')
save_2nifti(mask_G3_L_3d,f'ttest_G3_{ROI}_L_significantDiff_brainstat{opt2}')

mask_L_4d = np.stack([mask_G1_L_3d, mask_G2_L_3d], axis=-1) # , mask_G3_L_3d
save_2nifti(mask_L_4d,f'ttest_{ROI}_L_significantDiff_brainstat{opt2}_4d')

# R-thalamus
save_2nifti(mask_G1_R_3d,f'ttest_G1_{ROI}_R_significantDiff_brainstat{opt2}')
save_2nifti(mask_G2_R_3d,f'ttest_G2_{ROI}_R_significantDiff_brainstat{opt2}')
save_2nifti(mask_G3_R_3d,f'ttest_G3_{ROI}_R_significantDiff_brainstat{opt2}')

mask_R_4d = np.stack([mask_G1_R_3d, mask_G2_R_3d, mask_G3_R_3d], axis=-1)
save_2nifti(mask_R_4d,f'ttest_{ROI}_R_significantDiff_brainstat{opt2}_4d')


# save T-value
def reconstruct_tvalue2_3D(tval,outname):
    mask_3d = np.zeros(shape=roi_atlas.shape)
    mask_3d = np.reshape(mask_3d,nVoxels)
    for a, p in enumerate(par_val):
        t_par = tval.T[a]
        ind = np.where(roi_data == p)[0]
        mask_3d[ind] = t_par
    mask_3d = np.reshape(mask_3d,roi_atlas.shape)
    L_niftiImg = nib.Nifti1Image(mask_3d, roi_atlas.affine)
    nib.save(L_niftiImg,os.path.join(f'{out_dir}/ttest_{outname}.nii.gz'))
    return mask_3d

print(opt2)
G1_L_mask_tval = reconstruct_tvalue2_3D(G1_L_ttest,f'G1_{ROI}_L_thal_pmap_T_value{opt2}')
G2_L_mask_tval = reconstruct_tvalue2_3D(G2_L_ttest,f'G2_{ROI}_L_thal_pmap_T_value{opt2}')
G3_L_mask_tval = reconstruct_tvalue2_3D(G3_L_ttest,f'G3_{ROI}_L_thal_pmap_T_value{opt2}')

G1_R_mask_tval = reconstruct_tvalue2_3D(G1_R_ttest,f'G1_{ROI}_R_thal_pmap_T_value{opt2}')
G2_R_mask_tval = reconstruct_tvalue2_3D(G2_R_ttest,f'G2_{ROI}_R_thal_pmap_T_value{opt2}')
G3_R_mask_tval = reconstruct_tvalue2_3D(G3_R_ttest,f'G3_{ROI}_R_thal_pmap_T_value{opt2}')

tval_L_4d = np.stack([G1_L_mask_tval, G2_L_mask_tval, G3_L_mask_tval], axis=-1)
save_2nifti(tval_L_4d,f'ttest_{ROI}_L_thal_pmap_T_value{opt2}_4d')

tval_R_4d = np.stack([G1_R_mask_tval, G2_R_mask_tval, G3_R_mask_tval], axis=-1)
save_2nifti(tval_L_4d,f'ttest_{ROI}_R_thal_pmap_T_value{opt2}_4d')



###############################
##### 5. Figure
def avg_dict(G_L_data_thal):
    if isinstance(G_L_data_thal, dict):
        G_L_data_array = np.array(list(G_L_data_thal.values()))
    elif isinstance(G_L_data_thal, (list, np.ndarray)):
        G_L_data_array = G_L_data_thal
    else:
        ValueError("Input data must be either a dictionary or a NumPy array")
    G_L_data_avg = np.mean(G_L_data_array, axis=0)
    return G_L_data_avg

# L-thalamus
G1_L_data1_avg = avg_dict(G1_L_data1_pmap_list)
G1_L_data2_avg = avg_dict(G1_L_data2_pmap_list)

# R-thalamus
G1_R_data1_avg = avg_dict(G1_R_data1_pmap_list)
G1_R_data2_avg = avg_dict(G1_R_data2_pmap_list)



import matplotlib.pyplot as plt
import seaborn as sns

def plt_densMap(G_L_data1, G_L_data2, g, figsize=(5,5), xlim=(0.08, 0.08), save_path=None):
    fig, axs = plt.subplots(figsize=(figsize[0], figsize[1]))
    sns.kdeplot(G_L_data1, label='NOR', shade=True, color='#408abf') #'NOR' 'ROI 0'
    sns.kdeplot(G_L_data2, label='SCZ', shade=True, color='#ff6200') #'SCZ' 'ROI 1'
    # Add labels and title
    plt.xlabel(g, fontsize = 10)
    plt.ylabel('Density', fontsize = 10)
    axs.spines['top'].set_color('none')
    axs.spines['right'].set_color('none')
    axs.spines['left'].set_color('none')
    axs.set_yticks([])
    axs.set_yticklabels([])
    axs.set_xticks(xlim)
    axs.set_xlim(xlim[0], xlim[1])
    axs.tick_params(axis='x', labelsize=8)
    plt.legend(fontsize=8)
    plt.show()
    if save_path:
        fig.savefig(save_path, dpi=100, bbox_inches='tight')
        

ext = '_long'
plt_densMap(G1_L_data1_avg, G1_L_data2_avg, 'G1: L-thal projection', figsize=(8,2), xlim=(-1800,2500), save_path=f'{out_dir}/Figure/L_thal_grpdiff_distribution_G1_avg{ext}.png') # -1800,2500
plt_densMap(G1_R_data1_avg, G1_R_data2_avg, 'G1: R-thal projection', figsize=(8,2), xlim=(-1800, 2500), save_path=f'{out_dir}/Figure/R_thal_grpdiff_distribution_G1_avg{ext}.png') # -1800, 2500

