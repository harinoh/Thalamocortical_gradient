#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 08:52:08 2024

@author: harin_oh
"""

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import sys
import os
import statsmodels.formula.api as smf
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

roi= 'WholeCTX'
p=1000
ref_group='SCZOHC'
opt=f'_{ref_group}_FWHM4_Parcel{p}_IndivMul'
sub_group1 = 'OHC' # NOR
sub_group2 = 'SCZ' # patient
home_dir = '/Volumes/CCNC/harin_oh/1_Thalamocortical'
gradient_path = f'{home_dir}/1_Gradient_result/2_Alignment/SCZOHC/{roi}{opt}'
file_name_group1 = f'{home_dir}/SCZOHC/Subject_list_{sub_group1}.txt'
file_name_group2 = f'{home_dir}/SCZOHC/Subject_list_{sub_group2}.txt'
out_dir = f'{home_dir}/1_Gradient_result/3_1_Stat_anal_vox/SCZOHC/{roi}{opt}'

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
        L_thal = nib.load(f'{gradient_path}/{sub}_Gordon_L_thal_222_aligned_mul.cmaps.nii.gz')
        R_thal = nib.load(f'{gradient_path}/{sub}_Gordon_R_thal_222_aligned_mul.cmaps.nii.gz')
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
        
    return G1_L_thal, G2_L_thal, G3_L_thal, G1_R_thal, G2_R_thal, G3_R_thal


#### Main ###

if not os.path.exists(out_dir):
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
G1_L_data2, G2_L_data2, G3_L_data2, G1_R_data2, G2_R_data2, G3_R_data2 = nii_read(subject_list_group2)      # group SCZ

# Loading demographic data
cols = ['Folder_name', 'Group','age','sex_num','DOI'] #, 'sex_num'
demo_grp1 = pd.read_excel(f'{home_dir}/SCZOHC_rsfMRI_Demo.xlsx', sheet_name='baseline_n129', usecols=cols, nrows=69, index_col='Folder_name')
demo_grp2 = pd.read_excel(f'{home_dir}/SCZOHC_rsfMRI_Demo.xlsx', sheet_name='baseline_n129', usecols=cols, skiprows=range(1,70), nrows=60, index_col='Folder_name')
demo = pd.concat([demo_grp1,demo_grp2])

# Load FDRMS file
def read_FDRMS(subject_list):
    fdrms = {}
    for i in subject_list:
        fdrms_dir = f'{home_dir}/SCZOHC/{i}/preproc/{i}_MoCorr_FDRMS_val_MO'
        with open(fdrms_dir, 'r') as file:
            fdrms_handle = file.read().rstrip()
            fdrms_array = np.array(fdrms_handle.split("\n"), dtype=float)
            fdrms[i] = fdrms_array
    return fdrms

fdrms_grp1 = read_FDRMS(subject_list_group1)
fdrms_grp2 = read_FDRMS(subject_list_group2)

def mean_FDRMS(subject_list, fdrms):
    mean_fdrms = np.zeros(len(subject_list))
    for idx, sub in enumerate(subject_list):
        sub_fdrms = fdrms[sub]
        mean_fdrms[idx] = np.mean(sub_fdrms)
    return mean_fdrms

mean_fdrms_grp1 = mean_FDRMS(subject_list_group1, fdrms_grp1)
mean_fdrms_grp2 = mean_FDRMS(subject_list_group2, fdrms_grp2)

concat_mean_fdrms = np.concatenate((mean_fdrms_grp1,mean_fdrms_grp2))

# Load thalamic functional atlas
L_thal_atlas=nib.load(f'/{home_dir}/ROI_atlas/Gordon_Parcels/Gordon_L_thal_222.nii.gz')
R_thal_atlas=nib.load(f'/{home_dir}/ROI_atlas/Gordon_Parcels/Gordon_R_thal_222.nii.gz')
L_thal = L_thal_atlas.get_fdata()
R_thal = R_thal_atlas.get_fdata()
L_thal_indx = np.where(L_thal !=0)
R_thal_indx = np.where(R_thal !=0)

# Extract thalamic region
def thal_extract(G_data, thal_indx, subject_list):
    G_data_thal={}
    G_data_thal_list = []
    for sub in subject_list:
        sub_thal = G_data[sub][thal_indx[0], thal_indx[1], thal_indx[2]]
        G_data_thal[sub] = sub_thal
        G_data_thal_list.append(sub_thal)
    return G_data_thal, G_data_thal_list

G1_L_data1_thal, G1_L_data1_thal_list = thal_extract(G1_L_data1, L_thal_indx, subject_list_group1)
G2_L_data1_thal, G2_L_data1_thal_list = thal_extract(G2_L_data1, L_thal_indx, subject_list_group1)
G3_L_data1_thal, G3_L_data1_thal_list = thal_extract(G3_L_data1, L_thal_indx, subject_list_group1)
G1_R_data1_thal, G1_R_data1_thal_list = thal_extract(G1_R_data1, R_thal_indx, subject_list_group1)
G2_R_data1_thal, G2_R_data1_thal_list = thal_extract(G2_R_data1, R_thal_indx, subject_list_group1)
G3_R_data1_thal, G3_R_data1_thal_list = thal_extract(G3_R_data1, R_thal_indx, subject_list_group1)


G1_L_data2_thal, G1_L_data2_thal_list = thal_extract(G1_L_data2, L_thal_indx, subject_list_group2)
G2_L_data2_thal, G2_L_data2_thal_list = thal_extract(G2_L_data2, L_thal_indx, subject_list_group2)
G3_L_data2_thal, G3_L_data2_thal_list = thal_extract(G3_L_data2, L_thal_indx, subject_list_group2)
G1_R_data2_thal, G1_R_data2_thal_list = thal_extract(G1_R_data2, R_thal_indx, subject_list_group2)
G2_R_data2_thal, G2_R_data2_thal_list = thal_extract(G2_R_data2, R_thal_indx, subject_list_group2)
G3_R_data2_thal, G3_R_data2_thal_list = thal_extract(G3_R_data2, R_thal_indx, subject_list_group2)


#### Normalized version
# def norm_dict(G_L_data_thal):
#     norm_dic = {key : normalize([value], axis=1)[0] for key, value in G_L_data_thal.items()}
#     return norm_dic

# norm_G1_L_data1_thal = norm_dict(G1_L_data1_thal)
# norm_G2_L_data1_thal = norm_dict(G2_L_data1_thal)
# norm_G3_L_data1_thal = norm_dict(G3_L_data1_thal)
# norm_G1_R_data1_thal = norm_dict(G1_R_data1_thal)
# norm_G2_R_data1_thal = norm_dict(G2_R_data1_thal)
# norm_G3_R_data1_thal = norm_dict(G3_R_data1_thal)

# norm_G1_L_data2_thal = norm_dict(G1_L_data2_thal)
# norm_G2_L_data2_thal = norm_dict(G2_L_data2_thal)
# norm_G3_L_data2_thal = norm_dict(G3_L_data2_thal)
# norm_G1_R_data2_thal = norm_dict(G1_R_data2_thal)
# norm_G2_R_data2_thal = norm_dict(G2_R_data2_thal)
# norm_G3_R_data2_thal = norm_dict(G3_R_data2_thal)

## mean gradient value (for correlation)
def bi_avg_dict(G_L_data_thal,G_R_data_thal, subject_list):
    G_bi_data_avg = []
    for i in subject_list:
        sub_L_data = G_L_data_thal[i]
        sub_R_data = G_R_data_thal[i]
        bi_concat = np.concatenate((sub_L_data, sub_R_data))
        bi_max = max(bi_concat)
        bi_min = min(bi_concat)
        sub_avg = bi_max - bi_min
        # sub_avg = np.mean(bi_concat)
        G_bi_data_avg.append(sub_avg)
    return G_bi_data_avg

G1_Bi_data1_thal_avg = bi_avg_dict(G1_L_data1_thal, G1_R_data1_thal, subject_list_group1)
G2_Bi_data1_thal_avg = bi_avg_dict(G2_L_data1_thal, G2_R_data1_thal, subject_list_group1)
G3_Bi_data1_thal_avg = bi_avg_dict(G3_L_data1_thal, G3_R_data1_thal, subject_list_group1)

G1_Bi_data2_thal_avg = bi_avg_dict(G1_L_data2_thal, G1_R_data2_thal, subject_list_group2)
G2_Bi_data2_thal_avg = bi_avg_dict(G2_L_data2_thal, G2_R_data2_thal, subject_list_group2)
G3_Bi_data2_thal_avg = bi_avg_dict(G3_L_data2_thal, G3_R_data2_thal, subject_list_group2)

# correlation to fdrms (w/ NOR)
G1_corr_coef_grp1_fdrms = np.corrcoef(mean_fdrms_grp1, G1_Bi_data1_thal_avg)[0,1]
G2_corr_coef_grp1_fdrms = np.corrcoef(mean_fdrms_grp1, G2_Bi_data1_thal_avg)[0,1]
G3_corr_coef_grp1_fdrms = np.corrcoef(mean_fdrms_grp1, G3_Bi_data1_thal_avg)[0,1]
print(G1_corr_coef_grp1_fdrms, G2_corr_coef_grp1_fdrms, G3_corr_coef_grp1_fdrms)

# correlation to age (w/NOR)
G1_corr_coef_grp1_age = np.corrcoef(demo_grp1['age'], G1_Bi_data1_thal_avg)[0,1]
G2_corr_coef_grp1_age = np.corrcoef(demo_grp1['age'], G2_Bi_data1_thal_avg)[0,1]
G3_corr_coef_grp1_age = np.corrcoef(demo_grp1['age'], G3_Bi_data1_thal_avg)[0,1]
print(G1_corr_coef_grp1_age, G2_corr_coef_grp1_age, G3_corr_coef_grp1_age)

# correlation to sex (w/NOR)
G1_corr_coef_grp1_sex = np.corrcoef(demo_grp1['sex_num'], G1_Bi_data1_thal_avg)[0,1]
G2_corr_coef_grp1_sex = np.corrcoef(demo_grp1['sex_num'], G2_Bi_data1_thal_avg)[0,1]
G3_corr_coef_grp1_sex = np.corrcoef(demo_grp1['sex_num'], G3_Bi_data1_thal_avg)[0,1]
print(G1_corr_coef_grp1_sex, G2_corr_coef_grp1_sex, G3_corr_coef_grp1_sex)

# correlation to reduced volume (in SCZ)
# [G1_corr_coef_grp2_PULvol, G1_p_val_grp2_PULvol] = stats.pearsonr(VBM_pul_data2,G1_Bi_data2_thal_avg)
# [G2_corr_coef_grp2_PULvol, G2_p_val_grp2_PULvol] = stats.pearsonr(VBM_pul_data2,G2_Bi_data2_thal_avg)
# [G1_corr_coef_grp2_MDPULvol, G1_p_val_grp2_MDPULvol] = stats.pearsonr(VBM_md_pul_data2,G1_Bi_data2_thal_avg)
# [G2_corr_coef_grp2_MDPULvol, G2_p_val_grp2_MDPULvol] = stats.pearsonr(VBM_md_pul_data2,G2_Bi_data2_thal_avg)

# [G1_corr_coef_all_PULvol, G1_p_val_all_PULvol] = stats.pearsonr(VBM_pul_all,(G1_Bi_data1_thal_avg + G1_Bi_data2_thal_avg))
# [G2_corr_coef_all_PULvol, G2_p_val_all_PULvol] = stats.pearsonr(VBM_pul_all,(G1_Bi_data1_thal_avg + G2_Bi_data2_thal_avg))
# [G1_corr_coef_all_MDPULvol, G1_p_val_all_MDPULvol] = stats.pearsonr(VBM_md_pul_all,(G1_Bi_data1_thal_avg + G1_Bi_data2_thal_avg))
# [G2_corr_coef_all_MDPULvol, G2_p_val_all_MDPULvol] = stats.pearsonr(VBM_md_pul_all,(G1_Bi_data1_thal_avg + G2_Bi_data2_thal_avg))


# correlation to doi (w/NOR)
[G1_corr_coef_grp2_doi, G1_p_val_grp2_doi]  = np.corrcoef(demo_grp2['DOI'], G1_Bi_data2_thal_avg)
print(G1_corr_coef_grp2_doi, G1_p_val_grp2_doi)



###############################
#### 1. Regression
def regress_glm(G_L_data1_thal,G_L_data2_thal):
    vox_resid = []
    G_L_data1_list = list(G_L_data1_thal.values())
    G_L_data2_list = list(G_L_data2_thal.values())
    G_L_data1_arr = np.array(G_L_data1_list)
    G_L_data2_arr = np.array(G_L_data2_list)
        
    concat_array = np.row_stack((G_L_data1_arr, G_L_data2_arr))
        
    # Covariate + 1 designated to X
    X = np.column_stack((np.ones_like(demo['age']), demo['age'])) # demo['Group'], concat_mean_fdrms, demo['sex_num']
    model = LinearRegression()
    model.fit(X, concat_array)
    
    resid = concat_array - model.predict(X)
    vox_resid = resid
    vox_resid_data1 = vox_resid[0:69]
    vox_resid_data2 = vox_resid[69:129]
    return vox_resid_data1, vox_resid_data2

G1_L_resid_data1, G1_L_resid_data2 = regress_glm(G1_L_data1_thal,G1_L_data2_thal)
G1_R_resid_data1, G1_R_resid_data2 = regress_glm(G1_R_data1_thal,G1_R_data2_thal)




###############################
#### 2. Statistical testing
# format them into dataframe
demo['Group'] =  demo['Group'].astype('category')
demo['sex_num'] =  demo['sex_num'].astype('category')
subject_list_all = subject_list_group1 + subject_list_group2
df_meanFDRMS = pd.DataFrame(concat_mean_fdrms, columns=['meanFDRMS'], index=subject_list_all)
demo_FDRMS = pd.concat([demo, df_meanFDRMS], axis=1)


# t-test version
def ttest_construct(G_L_data1_thal, G_L_data2_thal, thal_ind, g):
    print(f'For {g}:')
    vox_t_stat = []
    vox_p_val = []
    for vox in range(0,len(thal_ind[0])):
        if isinstance(G_L_data1_thal, dict):
            G_L_data1_arr = np.array([v[vox] for k, v in G_L_data1_thal.items()])
            G_L_data2_arr = np.array([v[vox] for k, v in G_L_data2_thal.items()])
        elif isinstance(G_L_data1_thal, np.ndarray):
            G_L_data1_arr = G_L_data1_thal
            G_L_data2_arr = G_L_data2_thal
        else:
            ValueError("Input data must be either a dictionary or a NumPy array")
        
        t_stat, p_val = stats.ttest_ind(G_L_data1_arr,G_L_data2_arr,equal_var=False)
        vox_t_stat = (vox,t_stat)[1]
        vox_p_val = (vox, p_val)[1]
    
    uncorr_pval_find = np.where(vox_p_val < 0.05)[0]
    FDR_corrected_p_val = multipletests(vox_p_val, method='fdr_bh', alpha = 0.05)[1]
    FDR_corr_p_find = np.where(FDR_corrected_p_val < 0.05)[0]
    print('uncorrected voxels:', len(uncorr_pval_find))
    print('FDR_corr voxels:', len(FDR_corr_p_find))
    #print('surfstat qvalues:', len(qvalues))
    
    return vox_t_stat, vox_p_val, uncorr_pval_find, FDR_corrected_p_val, FDR_corr_p_find

# w/o regression
G1_L_t_stat, G1_L_p_val, G1_L_uncorr_p_find, G1_L_FDR_corr_p_val, G1_L_FDR_corr_p_find = ttest_construct(G1_L_data1_thal, G1_L_data2_thal, L_thal_indx, 'G1 L_thal')
G1_R_t_stat, G1_R_p_val, G1_R_uncorr_p_find, G1_R_FDR_corr_p_val, G1_R_FDR_corr_p_find = ttest_construct(G1_R_data1_thal, G1_R_data2_thal, R_thal_indx, 'G1 R_thal')

# w/ regression
G1_L_t_stat, G1_L_p_val, G1_L_uncorr_p_find, G1_L_FDR_corr_p_val_resid, G1_L_FDR_corr_p_find_resid = ttest_construct(G1_L_resid_data1, G1_L_resid_data2, L_thal_indx, 'G1 L_thal')

G1_R_t_stat, G1_R_p_val, G1_R_uncorr_p_find, G1_R_FDR_corr_p_val_resid, G1_R_FDR_corr_p_find_resid = ttest_construct(G1_R_resid_data1, G1_R_resid_data2, R_thal_indx, 'G1 R_thal')


## ANCOVA version
def ANCOVA_construct(G_L_data1_thal,G_L_data2_thal, roi, thal_ind, demo, G):
    print(f'Statistical result when {roi} ROI used')
    G_L_data1_df = pd.DataFrame(G_L_data1_thal).T
    G_L_data2_df = pd.DataFrame(G_L_data2_thal).T
    G_L_concat = pd.concat([G_L_data1_df, G_L_data2_df])
    f_stat = []
    p_val = []
    for vox in range(0,len(thal_ind[0])):
        G_vox = G_L_concat[[vox]].rename(columns={ vox : 'outcome'})
        G_concat = pd.concat([demo, G_vox],axis=1)
        # ANCOVA model fitting
        model = smf.ols('outcome ~ Group + age', data=G_concat).fit() # + meanFDRMS + sex_num
        f_stat.append(model.fvalue)
        p_val.append(model.f_pvalue)
    
    print('For', G)
    uncorr_p_find = np.where(np.array(p_val) < 0.05)[0]
    print('uncorrected voxels:', len(uncorr_p_find))
    Bonf_corr_p_val = multipletests(p_val, method='bonferroni', alpha = 0.05)[1]
    Bonf_corr_p_find = np.where(Bonf_corr_p_val < 0.05)[0]
    print('Bonf_corr voxels:' , len(Bonf_corr_p_find))
    
    FDR_corr_p_val = multipletests(p_val, method='fdr_bh', alpha = 0.05)[1]
    FDR_corr_p_find = np.where(FDR_corr_p_val < 0.05)[0]
    print('FDR_corr voxels:', len(FDR_corr_p_find))
    
    return f_stat, p_val, uncorr_p_find, FDR_corr_p_find#, model



G1_L_f_stat, G1_L_p_val, G1_L_uncorr_p_find, G1_L_FDR_corr_p_find = ANCOVA_construct(G1_L_data1_thal, G1_L_data2_thal, roi, L_thal_indx, demo_FDRMS, 'G1') # , G1_ancova_model

G1_R_f_stat, G1_R_p_val, G1_R_uncorr_p_find, G1_R_FDR_corr_p_find = ANCOVA_construct(G1_R_data1_thal, G1_R_data2_thal, roi, R_thal_indx, demo_FDRMS, 'G1')



## Surfstat version
from brainstat.stats.SLM import SLM
from brainstat.stats import terms

def brainstat_ttest(demo, G_L_data1_thal_list, G_L_data2_thal_list):
    term_intercept = terms.FixedEffect(1, names='intercept')
    term_age = terms.FixedEffect(demo['age'], names='age')
    term_group = terms.FixedEffect(demo['Group'], names='group')
    term_sex = terms.FixedEffect(demo['sex_num'], names='sex')
    #term_fdrms = terms.FixedEffect(concat_mean_fdrms)
    
    group = demo['Group']
    model = term_intercept + term_group + term_age #+ term_sex # + term_age + term_sex + term_fdrms 
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
    
    print('uncorrected p-val:', len(np.where(pval < 0.05)[0]))
    print('FDR corrected p-value:',len(np.where(qval < 0.05)[0]))
    
    del slm_age
    
    return ttest, qval, fdr_pval_idx, pval

G1_L_ttest, G1_L_qval, G1_L_qval_idx, G1_L_pval = brainstat_ttest(demo, G1_L_data1_thal_list, G1_L_data2_thal_list)
G1_R_ttest, G1_R_qval, G1_R_qval_idx, G1_R_pval = brainstat_ttest(demo, G1_R_data1_thal_list, G1_R_data2_thal_list)


##################################
## rebuild into 3D image to save
L_thal_affine = L_thal_atlas.affine
R_thal_affine = R_thal_atlas.affine

def save_nifti(thal_indx,FDR_corr_p_find,thal,thal_affine,outname, extra_FDR_corr_p_find=None):
    mask_thal = np.zeros(len(thal_indx[0]))
    mask_thal[FDR_corr_p_find] = 1
    mask_thal_3d = np.zeros(thal.shape)

    x = thal_indx[0]
    y = thal_indx[1]
    z = thal_indx[2]
    for k in range(0, len(thal_indx[0])):
        x_coord = x[k]
        y_coord = y[k]
        z_coord = z[k]
        mask_val = mask_thal[k]
        mask_thal_3d[x_coord, y_coord, z_coord] = mask_val
    
    nifti_img_3d = np.expand_dims(mask_thal_3d, axis=-1)
    
    if extra_FDR_corr_p_find is not None:
        extra_mask_thal = np.zeros(len(thal_indx[0]))
        extra_mask_thal[extra_FDR_corr_p_find] = 1
        extra_mask_thal_3d = np.zeros(thal.shape)
        
        x = thal_indx[0]
        y = thal_indx[1]
        z = thal_indx[2]
        for k in range(0, len(thal_indx[0])):
            x_coord = x[k]
            y_coord = y[k]
            z_coord = z[k]
            extra_mask_val = extra_mask_thal[k]
            extra_mask_thal_3d[x_coord, y_coord, z_coord] = extra_mask_val
        extra_nifti_img_3d = np.expand_dims(extra_mask_thal_3d, axis=-1)  # Add a new axis for the 4D image

        # Concatenate the 3D images along the new axis to create a 4D image
        nifti_img_4d = np.concatenate((nifti_img_3d, extra_nifti_img_3d), axis=-1)
    else:
        nifti_img_4d = nifti_img_3d
    
    nifti_img = nib.Nifti1Image(nifti_img_4d, thal_affine)
    nib.save(nifti_img,os.path.join(f'{out_dir}/ttest_{outname}.nii.gz'))
    
    return mask_thal_3d

G1_L_pval_idx=np.where(G1_L_pval<0.05)[0]
naming = 'brainstat_qvalue'
mask_G1_L_thal_3d = save_nifti(L_thal_indx,G1_L_qval_idx,L_thal,L_thal_affine,f'{naming}_G1_{roi}_L_thal_agesexCo')
mask_G1_R_thal_3d = save_nifti(R_thal_indx,G1_R_qval_idx,R_thal,R_thal_affine,f'{naming}_G1_{roi}_R_thal_agesexCo')

# uncorrected
mask_L_uncorrected = save_nifti(L_thal_indx, G1_L_uncorr_p_find, L_thal, L_thal_affine, f'{roi2}_L_thal_uncorr', G2_L_uncorr_p_find)

# Save T-value 
def save_T_val_2nifti(thal_indx,ttest,thal,outname):
    mask_thal = np.zeros(len(thal_indx[0]))
    mask_thal = ttest
    mask_thal_3d = np.zeros(thal.shape)
    
    x = thal_indx[0]
    y = thal_indx[1]
    z = thal_indx[2]
    for k in range(0, len(thal_indx[0])):
        x_coord = x[k]
        y_coord = y[k]
        z_coord = z[k]
        mask_val = mask_thal[0][k]
        mask_thal_3d[x_coord, y_coord, z_coord] = mask_val
    L_niftiImg = nib.Nifti1Image(mask_thal_3d, L_thal_affine)
    nib.save(L_niftiImg,os.path.join(f'{out_dir}/ttest_{outname}.nii.gz'))

    return mask_thal_3d

G1_L_thal_Tval = save_T_val_2nifti(L_thal_indx, G1_L_ttest, L_thal, f'G1_{roi}_L_thal_T_value')
G2_L_thal_Tval = save_T_val_2nifti(L_thal_indx, G2_L_ttest, L_thal, f'G2_{roi}_L_thal_T_value')

G1_R_thal_Tval = save_T_val_2nifti(R_thal_indx, G1_R_ttest, R_thal, f'G1_{roi}_R_thal_T_value')
G2_R_thal_Tval = save_T_val_2nifti(R_thal_indx, G2_R_ttest, R_thal, f'G2_{roi}_R_thal_T_value')

G1_Bi_thal_Tval = G1_L_thal_Tval + G1_R_thal_Tval
G2_Bi_thal_Tval = G2_L_thal_Tval + G2_R_thal_Tval

G1_thalImg = nib.Nifti1Image(G1_Bi_thal_Tval, L_thal_affine)
nib.save(G1_thalImg,os.path.join(f'{out_dir}/ttest_G1_{roi}_Bi_thal_T_value.nii.gz'))
G2_thalImg = nib.Nifti1Image(G2_Bi_thal_Tval, L_thal_affine)
nib.save(G2_thalImg,os.path.join(f'{out_dir}/ttest_G2_{roi}_Bi_thal_T_value.nii.gz'))





###############################
#### Distribution graph
# average ver
def avg_dict(G_L_data_thal):
    if isinstance(G_L_data_thal, dict):
        G_L_data_array = np.array(list(G_L_data_thal.values()))
    elif isinstance(G_L_data_thal, (list, np.ndarray)):
        G_L_data_array = G_L_data_thal
    else:
        ValueError("Input data must be either a dictionary or a NumPy array")
    G_L_data_avg = np.mean(G_L_data_array, axis=0)
    return G_L_data_avg

G1_L_data1_avg = avg_dict(G1_L_data1_thal)
G1_R_data1_avg = avg_dict(G1_R_data1_thal)

G1_L_data2_avg = avg_dict(G1_L_data2_thal)
G1_R_data2_avg = avg_dict(G1_R_data2_thal)

# regress version
G1_L_data1_avg = avg_dict(G1_L_resid_data1)
G1_R_data1_avg = avg_dict(G1_R_resid_data1)

G1_L_data2_avg = avg_dict(G1_L_resid_data2)
G1_R_data2_avg = avg_dict(G1_R_resid_data2)

# #norm
# G1_L_data1_avg_norm = avg_dict(norm_G1_L_data1_thal)
# G2_L_data1_avg_norm = avg_dict(norm_G2_L_data1_thal)

# G1_L_data2_avg_norm = avg_dict(norm_G1_L_data2_thal)
# G2_L_data2_avg_norm = avg_dict(norm_G2_L_data2_thal)

# flatten ver
def flatten_data(G_L_data_thal):
    G_data_flatten = []
    if isinstance(G_L_data_thal, dict):
        for key in G_L_data_thal:
            G_data_flatten.extend(G_L_data_thal[key])
    elif isinstance(G_L_data_thal, np.ndarray):
        G_data_flatten = G_L_data_thal.reshape(-1)
    return G_data_flatten

G1_L_data1_flatten = flatten_data(G1_L_data1_thal)
G2_L_data1_flatten = flatten_data(G2_L_data1_thal)
G3_L_data1_flatten = flatten_data(G3_L_data1_thal)
G1_R_data1_flatten = flatten_data(G1_R_data1_thal)
G2_R_data1_flatten = flatten_data(G2_R_data1_thal)
G3_R_data1_flatten = flatten_data(G3_R_data1_thal)

G1_L_data2_flatten = flatten_data(G1_L_data2_thal)
G2_L_data2_flatten = flatten_data(G2_L_data2_thal)
G3_L_data2_flatten = flatten_data(G3_L_data2_thal)
G1_R_data2_flatten = flatten_data(G1_R_data2_thal)
G2_R_data2_flatten = flatten_data(G2_R_data2_thal)
G3_R_data2_flatten = flatten_data(G3_R_data2_thal)

# regress version
G1_L_data1_flatten = flatten_data(G1_L_resid_data1)
G2_L_data1_flatten = flatten_data(G2_L_resid_data1)
G3_L_data1_flatten = flatten_data(G3_L_resid_data1)
G1_R_data1_flatten = flatten_data(G1_R_resid_data1)
G2_R_data1_flatten = flatten_data(G2_R_resid_data1)
G3_R_data1_flatten = flatten_data(G3_R_resid_data1)

G1_L_data2_flatten = flatten_data(G1_L_resid_data2)
G2_L_data2_flatten = flatten_data(G2_L_resid_data2)
G3_L_data2_flatten = flatten_data(G3_L_resid_data2)
G1_R_data2_flatten = flatten_data(G1_R_resid_data2)
G2_R_data2_flatten = flatten_data(G2_R_resid_data2)
G3_R_data2_flatten = flatten_data(G3_R_resid_data2)

# concatenate L & R thalamus
G1_Bi_data1_flat = G1_L_data1_flatten + G1_R_data1_flatten
G2_Bi_data1_flat = G2_L_data1_flatten + G2_R_data1_flatten
G3_Bi_data1_flat = G3_L_data1_flatten + G3_R_data1_flatten

G1_Bi_data2_flat = G1_L_data2_flatten + G1_R_data2_flatten
G2_Bi_data2_flat = G2_L_data2_flatten + G2_R_data2_flatten
G3_Bi_data2_flat = G3_L_data2_flatten + G3_R_data2_flatten


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

ext = '_regress_sexCo'
ext = '_long'
plt_densMap(G1_Bi_data1_flat, G1_Bi_data2_flat, 'G1: Bi-thal', figsize=(8,4), xlim=(-0.0085, 0.006), save_path=f'{out_dir}/Figure/Bi_thal_grpdiff_distribution_G1_flat{ext}.png')

plt_densMap(G1_L_data1_flatten, G1_L_data2_flatten, 'G1: L-thal', figsize=(8,4), xlim=(-0.0085, 0.006), save_path=f'{out_dir}/Figure/L_thal_grpdiff_distribution_G1_flat{ext}.png')





def dict_concat(G_L_data_thal, G_R_data_thal):
    G1_Bi_data = {}
    for key in G_L_data_thal.keys():
        G1_Bi_data[key] = np.concatenate([G_L_data_thal[key], G_R_data_thal[key]])
        return G1_Bi_data

G1_Bi_data1 = dict_concat(G1_L_data1_thal, G1_R_data1_thal) # concatenate L & R thalamus
G1_Bi_data2 = dict_concat(G1_L_data2_thal, G1_R_data2_thal)

G1_Bi_data1_avg = avg_dict(G1_Bi_data1)
G1_Bi_data2_avg = avg_dict(G1_Bi_data2)
