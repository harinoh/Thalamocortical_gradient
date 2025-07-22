#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 20:47:43 2024

@author: harin_oh
"""

import nibabel as nib
import os
import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.linear_model import LinearRegression
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
from statsmodels.stats.multitest import multipletests

group = 'SCZOHC'
roi= 'WholeCTX'
ref_group='SCZOHC'
opt='_FWHM4_noRearrange_7gradient'
sub_group1 = 'OHC' # NOR
sub_group2 = 'SCZ' # patient
gradient_path = f'/Volume/CCNC/harin_oh/1_thalamocortical/Gradient_result/2_Alignment/{group}/{roi}_{ref_group}{opt}'
file_name_group1 = f'/Volume/CCNC/harin_oh/1_thalamocortical/{group}/Subject_list_{sub_group1}.txt'
file_name_group2 = f'/Volume/CCNC/harin_oh/1_thalamocortical/{group}/Subject_list_{sub_group2}.txt'
home_dir = '/Volume/CCNC/harin_oh/1_thalamocortical'
out_dir = f'{home_dir}/Gradient_result/3_4_SymptomCorr/{group}/{roi}_{ref_group}{opt}'

def read_subject_list(file_name):
    with open(file_name, 'r') as my_file:
        file_handle = my_file.read()
        subjects = file_handle.rstrip()
        subject_list = subjects.split("\n")
    return subject_list

def nii_read(subject_list):
    G1 = {}
    G2 = {}
    G3 = {}
    for i in subject_list:
        L_thal = nib.load(f'{gradient_path}/{i}_Gordon_L_thal_222_aligned.cmaps.nii.gz')
        R_thal = nib.load(f'{gradient_path}/{i}_Gordon_R_thal_222_aligned.cmaps.nii.gz')
        Bi_thal_G1 = L_thal.get_fdata()[:, :, :, 0] + R_thal.get_fdata()[:, :, :, 0]
        Bi_thal_G2 = L_thal.get_fdata()[:, :, :, 1] + R_thal.get_fdata()[:, :, :, 1]
        Bi_thal_G3 = L_thal.get_fdata()[:, :, :, 2] + R_thal.get_fdata()[:, :, :, 2]
        
        G1[i] = Bi_thal_G1
        G2[i] = Bi_thal_G2
        G3[i] = Bi_thal_G3
        
    return G1, G2, G3

def rem_zero_val(subject_list, data, Thal_ind, Thal_x, Thal_y, Thal_z):
    G_nonzero_data = {}    
    # Remove zero values
    for sub in subject_list:
        G_nonzero_data[sub] = np.zeros((1, len(Thal_ind)))
        for k in range(0, len(Thal_ind)):
            thal_x_coord = thal_x[k]
            thal_y_coord = thal_y[k]
            thal_z_coord = thal_z[k]
            nonzero_val = data[sub][thal_x_coord, thal_y_coord, thal_z_coord]
            G_nonzero_data[sub][:,k] = nonzero_val
    return G_nonzero_data

def regress_out_age_sex(df_demo, G_list):
    X = np.column_stack((np.ones_like(df_demo['age']), df_demo['age']))
    model = LinearRegression().fit(X, G_list)
    predicted = model.predict(X)
    G_adjust = G_list - predicted
    return G_adjust

if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    print(f'Folder {out_dir} created.')
else:
    print(f'from sklearn import preprocessingder {out_dir} already exists')


#### 0. Load data
# read demographic data
cols = ['Hospital ID', 'Folder_name', 'age', 'sex_num', 'IQ', 'B_PANSS_T', 'B_PANSS_P', 'B_PANSS_N', 'B_PANSS_G', 'TMT_A', 'TMT_B', 'Verbal_correct', 'Visual_correct', 'LNS', 'Spatial Span']
df_demo = pd.read_excel(f'{home_dir}/BCS_ECT_OHC_rsfMRI_Demo.xlsx', sheet_name='baseline_n129', usecols=cols)
df_grp1_demo = df_demo.iloc[0:69]
df_grp2_demo = df_demo.iloc[69:129]

# read subject list file
subject_list_group1 = read_subject_list(file_name_group1)
subject_list_group2 = read_subject_list(file_name_group2)

G1_data1, G2_data1, G3_data1 = nii_read(subject_list_group1)      # group NOR
G1_data2, G2_data2, G3_data2 = nii_read(subject_list_group2)      # group SCZ

# load thalamic template
Gordon_Thal = nib.load(f'{home_dir}/ROI_atlas/Gordon_Parcels/Gordon_Bi_thal_222.nii.gz').get_fdata()
Gordon_thal_coord = Gordon_Thal !=0
Gordon_Thal_ind = np.argwhere(Gordon_thal_coord)
thal_x, thal_y, thal_z = Gordon_Thal_ind[:,0], Gordon_Thal_ind[:,1], Gordon_Thal_ind[:,2]

# extract thalamic area
G1_nonzero_data1 = rem_zero_val(subject_list_group1,G1_data1, Gordon_Thal_ind, thal_x, thal_y, thal_z)
G1_nonzero_data2 = rem_zero_val(subject_list_group2,G1_data2, Gordon_Thal_ind, thal_x, thal_y, thal_z)

# convert to list
G1_list_data1 = np.array([v[0] for k, v in G1_nonzero_data1.items()])
G1_list_data2 = np.array([v[0] for k, v in G1_nonzero_data2.items()])



###############################
#### 1. prepare data (regression and/or normalization)
# Regressing out covariates
G1_adjust_data1 = regress_out_age_sex(df_grp1_demo, G1_list_data1)
G1_adjust_data2 = regress_out_age_sex(df_grp2_demo, G1_list_data2)

# getting the mean gradient values of each individual - unnormalized ver.
G1_mean_data1 = np.mean(G1_adjust_data1,axis=1)
G1_mean_data2 = np.mean(G1_adjust_data2,axis=1)


# w/o regression mean gradient
G1_mean_data1 = np.mean(G1_list_data1,axis=1)
G1_mean_data2 = np.mean(G1_list_data2,axis=1)


# normalize data - w/o regression
G1_norm_data1 = preprocessing.normalize(G1_list_data1)
G1_norm_data2 = preprocessing.normalize(G1_list_data2)


# getting the mean gradient values of each individual - normalized ver.
G1_mean_data1 = np.mean(G1_norm_data1, axis=1)
G1_mean_data2 = np.mean(G1_norm_data2, axis=1)




###############################
#### 2. Correlation to symptoms 
## symptom scores only compared in patient group
# PANSS_total
df_grp2_PANSS_T = df_grp2_demo['B_PANSS_T'].tolist()
[corr_coeff_G1_panssT, corr_pval_G1_panssT] = stats.pearsonr(G1_mean_data2, df_grp2_PANSS_T)
[corr_coeff_G2_panssT, corr_pval_G2_panssT] = stats.pearsonr(G2_mean_data2, df_grp2_PANSS_T)

# PANSS_positive
df_grp2_PANSS_P = df_grp2_demo['B_PANSS_P'].tolist()
[corr_coeff_G1_panssP, corr_pval_G1_panssP] = stats.pearsonr(G1_mean_data2, df_grp2_PANSS_P)
[corr_coeff_G2_panssP, corr_pval_G2_panssP] = stats.pearsonr(G2_mean_data2, df_grp2_PANSS_P)

# PANSS_negative
df_grp2_PANSS_N = df_grp2_demo['B_PANSS_N'].tolist()
[corr_coeff_G1_panssN, corr_pval_G1_panssN] = stats.pearsonr(G1_mean_data2, df_grp2_PANSS_N)
[corr_coeff_G2_panssN, corr_pval_G2_panssN] = stats.pearsonr(G2_mean_data2, df_grp2_PANSS_N)

# PANSS_general
df_grp2_PANSS_G = df_grp2_demo['B_PANSS_G'].tolist()
[corr_coeff_G1_panssG, corr_pval_G1_panssG] = stats.pearsonr(G1_mean_data2, df_grp2_PANSS_G)
[corr_coeff_G2_panssG, corr_pval_G2_panssG] = stats.pearsonr(G2_mean_data2, df_grp2_PANSS_G)

# concat all p_val & multiple correction
PANSS_pval_G1 = np.array([corr_pval_G1_panssT, corr_pval_G1_panssP, corr_pval_G1_panssN, corr_pval_G1_panssG])
PANSS_pval_G2 = np.array([corr_pval_G2_panssT, corr_pval_G2_panssP, corr_pval_G2_panssN, corr_pval_G2_panssG])
FDR_corrected_p_val_G1 = multipletests(PANSS_pval_G1, method='fdr_bh', alpha = 0.05)[1]
FDR_corrected_p_val_G2 = multipletests(PANSS_pval_G2, method='fdr_bh', alpha = 0.05)[1]

print('Gradient & symptom association - normalized version')
print('         correlation coefficient    p-value    corr_p-val')
print(f'PANSS_T: G1 ','{:.5f}'.format(corr_coeff_G1_panssT),'    ','{:.5f}'.format(corr_pval_G1_panssT),'    ','{:.5f}'.format(FDR_corrected_p_val_G1[0])) 
print(f'         G2 ','{:.5f}'.format(corr_coeff_G2_panssT),'    ','{:.5f}'.format(corr_pval_G2_panssT),'    ','{:.5f}'.format(FDR_corrected_p_val_G2[0])) 
print(f'PANSS_P: G1 ','{:.5f}'.format(corr_coeff_G1_panssP),'    ','{:.5f}'.format(corr_pval_G1_panssP),'    ','{:.5f}'.format(FDR_corrected_p_val_G1[1])) 
print(f'         G2 ','{:.5f}'.format(corr_coeff_G2_panssP),'    ','{:.5f}'.format(corr_pval_G2_panssP),'    ','{:.5f}'.format(FDR_corrected_p_val_G2[1])) 
print(f'PANSS_N: G1 ','{:.5f}'.format(corr_coeff_G1_panssN),'    ','{:.5f}'.format(corr_pval_G1_panssN),'    ','{:.5f}'.format(FDR_corrected_p_val_G1[2])) 
print(f'         G2 ','{:.5f}'.format(corr_coeff_G2_panssN),'    ','{:.5f}'.format(corr_pval_G2_panssN),'    ','{:.5f}'.format(FDR_corrected_p_val_G2[2])) 
print(f'PANSS_G: G1 ','{:.5f}'.format(corr_coeff_G1_panssG),'    ','{:.5f}'.format(corr_pval_G1_panssG),'    ','{:.5f}'.format(FDR_corrected_p_val_G1[3])) 
print(f'         G2 ','{:.5f}'.format(corr_coeff_G2_panssG),'    ','{:.5f}'.format(corr_pval_G2_panssG),'    ','{:.5f}'.format(FDR_corrected_p_val_G2[3])) 

#_normalized version
G1_norm_grp1_T = np.transpose(G1_norm_data1)
G1_norm_grp2_T = np.transpose(G1_norm_data2)
[corr_coeff_G1_panssT, corr_pval_G1_panssT] = stats.pearsonr(G1_norm_grp1_T, G1_norm_grp2_T)




###############################
#### 3. Correlation to Cognition
def cog_corr(df_grp_demo, cog_column, G1_mean_data, G2_mean_data):
    df_grp_cog = np.array(df_grp_demo[cog_column])
    grp_cog_null_indx = np.where(~np.isnan(df_grp_cog))[0] # Get indices as an array
    grp_cog_sub = df_grp_cog[grp_cog_null_indx]
    G1_mean_data_sub = G1_mean_data[grp_cog_null_indx]
    G2_mean_data_sub = G2_mean_data[grp_cog_null_indx]
    n = len(grp_cog_null_indx)
    [corr_coeff_G1_cog_grp, corr_pval_G1_cog_grp] = stats.pearsonr(G1_mean_data_sub, grp_cog_sub)
    [corr_coeff_G2_cog_grp, corr_pval_G2_cog_grp] = stats.pearsonr(G2_mean_data_sub, grp_cog_sub)
    return corr_coeff_G1_cog_grp, corr_pval_G1_cog_grp, corr_coeff_G2_cog_grp, corr_pval_G2_cog_grp, n, grp_cog_sub

## IQ 
corr_coeff_G1_IQ_grp1, corr_pval_G1_IQ_grp1, corr_coeff_G2_IQ_grp1, corr_pval_G2_IQ_grp1, n_IQ_grp1, grp1_IQ_demo = cog_corr(df_grp1_demo, 'IQ', G1_mean_data1, G2_mean_data1)
corr_coeff_G1_IQ_grp2, corr_pval_G1_IQ_grp2, corr_coeff_G2_IQ_grp2, corr_pval_G2_IQ_grp2, n_IQ_grp2 = cog_corr(df_grp2_demo, 'IQ', G1_mean_data2, G2_mean_data2)

## TMT 
corr_coeff_G1_TMT_A_grp1, corr_pval_G1_TMT_A_grp1, corr_coeff_G2_TMT_A_grp1, corr_pval_G2_TMT_A_grp1, n_TMT_A_grp1, grp1_TMT_A_demo = cog_corr(df_grp1_demo, 'TMT_A', G1_mean_data1, G2_mean_data1)
corr_coeff_G1_TMT_A_grp2, corr_pval_G1_TMT_A_grp2, corr_coeff_G2_TMT_A_grp2, corr_pval_G2_TMT_A_grp2, n_TMT_A_grp2, grp2_TMT_A_demo = cog_corr(df_grp2_demo, 'TMT_A', G1_mean_data2, G2_mean_data2)

corr_coeff_G1_TMT_B_grp1, corr_pval_G1_TMT_B_grp1, corr_coeff_G2_TMT_B_grp1, corr_pval_G2_TMT_B_grp1, n_TMT_B_grp1, grp1_TMT_B_demo = cog_corr(df_grp1_demo, 'TMT_B', G1_mean_data1, G2_mean_data1)
corr_coeff_G1_TMT_B_grp2, corr_pval_G1_TMT_B_grp2, corr_coeff_G2_TMT_B_grp2, corr_pval_G2_TMT_B_grp2, n_TMT_B_grp2, grp2_TMT_B_demo = cog_corr(df_grp2_demo, 'TMT_B', G1_mean_data2, G2_mean_data2)

# Verbal/ Visual FT 
corr_coeff_G1_verbalFT_grp1, corr_pval_G1_verbalFT_grp1, corr_coeff_G2_verbalFT_grp1, corr_pval_G2_verbalFT_grp1, n_verbalFT_grp1, grp1_VerbalFT_demo = cog_corr(df_grp1_demo, 'Verbal_correct', G1_mean_data1, G2_mean_data1)
corr_coeff_G1_verbalFT_grp2, corr_pval_G1_verbalFT_grp2, corr_coeff_G2_verbalFT_grp2, corr_pval_G2_verbalFT_grp2, n_verbalFT_grp2, grp2_VerbalFT_demo = cog_corr(df_grp2_demo, 'Verbal_correct', G1_mean_data2, G2_mean_data2)

corr_coeff_G1_visualFT_grp1, corr_pval_G1_visualFT_grp1, corr_coeff_G2_visualFT_grp1, corr_pval_G2_visualFT_grp1, n_visualFT_grp1, grp1_VisualFT_demo = cog_corr(df_grp1_demo, 'Visual_correct', G1_mean_data1, G2_mean_data1)
corr_coeff_G1_visualFT_grp2, corr_pval_G1_visualFT_grp2, corr_coeff_G2_visualFT_grp2, corr_pval_G2_visualFT_grp2, n_visualFT_grp2, grp1_VisualFT_demo = cog_corr(df_grp2_demo, 'Visual_correct', G1_mean_data2, G2_mean_data2)

# LNS/Corsi
corr_coeff_G1_LNS_grp1, corr_pval_G1_LNS_grp1, corr_coeff_G2_LNS_grp1, corr_pval_G2_LNS_grp1, n_LNS_grp1, grp1_LNS_demo = cog_corr(df_grp1_demo, 'LNS', G1_mean_data1, G2_mean_data1)
corr_coeff_G1_LNS_grp2, corr_pval_G1_LNS_grp2, corr_coeff_G2_LNS_grp2, corr_pval_G2_LNS_grp2, n_LNS_grp2, grp2_LNS_demo = cog_corr(df_grp2_demo, 'LNS', G1_mean_data2, G2_mean_data2)

corr_coeff_G1_Corsi_grp1, corr_pval_G1_Corsi_grp1, corr_coeff_G2_Corsi_grp1, corr_pval_G2_Corsi_grp1, n_Corsi_grp1, grp1_Corsi_demo = cog_corr(df_grp1_demo, 'Spatial Span', G1_mean_data1, G2_mean_data1)
corr_coeff_G1_Corsi_grp2, corr_pval_G1_Corsi_grp2, corr_coeff_G2_Corsi_grp2, corr_pval_G2_Corsi_grp2, n_Corsi_grp2, grp2_Corsi_demo = cog_corr(df_grp2_demo, 'Spatial Span', G1_mean_data2, G2_mean_data2)

# Multiple test correction
# NOR
G1_grp1_pval = np.array([corr_pval_G1_IQ_grp1, corr_pval_G1_TMT_A_grp1, corr_pval_G1_TMT_B_grp1, corr_pval_G1_verbalFT_grp1,  corr_pval_G1_LNS_grp1, corr_pval_G1_Corsi_grp1], dtype=np.float64) #corr_pval_G1_visualFT_grp1,
G2_grp1_pval = np.array([corr_pval_G2_IQ_grp1, corr_pval_G2_TMT_A_grp1, corr_pval_G2_TMT_B_grp1, corr_pval_G2_verbalFT_grp1,  corr_pval_G2_LNS_grp1, corr_pval_G2_Corsi_grp1], dtype=np.float64) #corr_pval_G2_visualFT_grp1,

G1_grp2_pval = np.array([corr_pval_G1_IQ_grp2, corr_pval_G1_TMT_A_grp2, corr_pval_G1_TMT_B_grp2, corr_pval_G1_verbalFT_grp2,  corr_pval_G1_LNS_grp2, corr_pval_G1_Corsi_grp2], dtype=np.float64) #corr_pval_G1_visualFT_grp2,
G2_grp2_pval = np.array([corr_pval_G2_IQ_grp2, corr_pval_G2_TMT_A_grp2, corr_pval_G2_TMT_B_grp2, corr_pval_G2_verbalFT_grp2, corr_pval_G2_LNS_grp2, corr_pval_G2_Corsi_grp2], dtype=np.float64) #corr_pval_G2_visualFT_grp2,

FDR_corr_G1_grp1 = multipletests(G1_grp1_pval, method='fdr_bh', alpha = 0.05)[1]
FDR_corr_G2_grp1 = multipletests(G2_grp1_pval, method='fdr_bh', alpha = 0.05)[1]

FDR_corr_G1_grp2 = multipletests(G1_grp2_pval, method='fdr_bh', alpha = 0.05)[1]
FDR_corr_G2_grp2 = multipletests(G2_grp2_pval, method='fdr_bh', alpha = 0.05)[1]


## Print result - def function analyze G1 & G2 of a single group
def print_corr_results(sub_group, cog_name, n_cog_grp, corr_coeff_G1_cog_grp, corr_pval_G1_cog_grp, FDR_corr_G1_grp, corr_coeff_G2_cog_grp, corr_pval_G2_cog_grp, FDR_corr_G2_grp, col_num):
    print(f'{sub_group} - {cog_name}:  G1 {n_cog_grp}   ', '{:.5f}'.format(corr_coeff_G1_cog_grp), '     ', '{:.5f}'.format(corr_pval_G1_cog_grp), '     ', '{:.5f}'.format(FDR_corr_G1_grp[col_num]))
    print(f'           G2     ', '{:.5f}'.format(corr_coeff_G2_cog_grp), '     ', '{:.5f}'.format(corr_pval_G2_cog_grp), '     ', '{:.5f}'.format(FDR_corr_G2_grp[col_num]))

print('             n      corr coeff    p-value      corr p-value')
print_corr_results(sub_group1, 'IQ', n_IQ_grp1, corr_coeff_G1_IQ_grp1, corr_pval_G1_IQ_grp1, FDR_corr_G1_grp1, corr_coeff_G2_IQ_grp1, corr_pval_G2_IQ_grp1, FDR_corr_G2_grp1, 0)
print_corr_results(sub_group2, 'IQ', n_IQ_grp2, corr_coeff_G1_IQ_grp2, corr_pval_G1_IQ_grp2, FDR_corr_G1_grp2, corr_coeff_G2_IQ_grp2, corr_pval_G2_IQ_grp2, FDR_corr_G2_grp2, 0)

print_corr_results(sub_group1, 'TMT_A', n_TMT_A_grp1, corr_coeff_G1_TMT_A_grp1, corr_pval_G1_TMT_A_grp1, FDR_corr_G1_grp1, corr_coeff_G2_TMT_A_grp1, corr_pval_G2_TMT_A_grp1, FDR_corr_G2_grp1, 1)
print_corr_results(sub_group2, 'TMT_A', n_TMT_A_grp2, corr_coeff_G1_TMT_A_grp2, corr_pval_G1_TMT_A_grp2, FDR_corr_G1_grp2, corr_coeff_G2_TMT_A_grp2, corr_pval_G2_TMT_A_grp2, FDR_corr_G2_grp2, 1)

print_corr_results(sub_group1, 'TMT_B', n_TMT_B_grp1, corr_coeff_G1_TMT_B_grp1, corr_pval_G1_TMT_B_grp1, FDR_corr_G1_grp1, corr_coeff_G2_TMT_B_grp1, corr_pval_G2_TMT_B_grp1, FDR_corr_G2_grp1, 2)
print_corr_results(sub_group2, 'TMT_B', n_TMT_B_grp2, corr_coeff_G1_TMT_B_grp2, corr_pval_G1_TMT_B_grp2, FDR_corr_G1_grp2, corr_coeff_G2_TMT_B_grp2, corr_pval_G2_TMT_B_grp2, FDR_corr_G2_grp2, 2)

print_corr_results(sub_group1, 'verbalFT', n_verbalFT_grp1, corr_coeff_G1_verbalFT_grp1, corr_pval_G1_verbalFT_grp1, FDR_corr_G1_grp1, corr_coeff_G2_verbalFT_grp1, corr_pval_G2_verbalFT_grp1, FDR_corr_G2_grp1, 3)
print_corr_results(sub_group2, 'verbalFT', n_verbalFT_grp2, corr_coeff_G1_verbalFT_grp2, corr_pval_G1_verbalFT_grp2, FDR_corr_G1_grp2, corr_coeff_G2_verbalFT_grp2, corr_pval_G2_verbalFT_grp2, FDR_corr_G2_grp2, 3)

print_corr_results(sub_group1, 'visualFT', n_visualFT_grp1, corr_coeff_G1_visualFT_grp1, corr_pval_G1_visualFT_grp1, FDR_corr_G1_grp1, corr_coeff_G2_visualFT_grp1, corr_pval_G2_visualFT_grp1, FDR_corr_G2_grp1, 4)
print_corr_results(sub_group2, 'visualFT', n_visualFT_grp2, corr_coeff_G1_visualFT_grp2, corr_pval_G1_visualFT_grp2, FDR_corr_G1_grp2, corr_coeff_G2_visualFT_grp2, corr_pval_G2_visualFT_grp2, FDR_corr_G2_grp2, 4)

print_corr_results(sub_group1, 'LNS', n_LNS_grp1, corr_coeff_G1_LNS_grp1, corr_pval_G1_LNS_grp1, FDR_corr_G1_grp1, corr_coeff_G2_LNS_grp1, corr_pval_G2_LNS_grp1, FDR_corr_G2_grp1, 5)
print_corr_results(sub_group2, 'LNS', n_LNS_grp2, corr_coeff_G1_LNS_grp2, corr_pval_G1_LNS_grp2, FDR_corr_G1_grp2, corr_coeff_G2_LNS_grp2, corr_pval_G2_LNS_grp2, FDR_corr_G2_grp2, 5)

print_corr_results(sub_group1, 'Corsi', n_Corsi_grp1, corr_coeff_G1_Corsi_grp1, corr_pval_G1_Corsi_grp1, FDR_corr_G1_grp1, corr_coeff_G2_Corsi_grp1, corr_pval_G2_Corsi_grp1, FDR_corr_G2_grp1, 6)
print_corr_results(sub_group2, 'Corsi', n_Corsi_grp2, corr_coeff_G1_Corsi_grp2, corr_pval_G1_Corsi_grp2, FDR_corr_G1_grp2, corr_coeff_G2_Corsi_grp2, corr_pval_G2_Corsi_grp2, FDR_corr_G2_grp2, 6)

## mean & std
import statistics

def print_mean_sd(cog,df_grp_demo):
    grp_cog_arr = np.array(df_grp_demo[cog])
    grp_cog_null_indx = np.where(~np.isnan(grp_cog_arr))[0] # Get indices as an array
    grp_cog_sub = grp_cog_arr[grp_cog_null_indx]
    grp_cog_mean = statistics.mean(grp_cog_sub)
    grp_cog_sd = statistics.stdev(grp_cog_sub)
    print(f'{cog}:   ','{:.1f}'.format(grp_cog_mean),'    ','{:.2f}'.format(grp_cog_sd))

print('NOR grp       mean        std')
print_mean_sd('IQ',df_grp1_demo)
print_mean_sd('TMT_A',df_grp1_demo)
print_mean_sd('TMT_B',df_grp1_demo)
print_mean_sd('Verbal_correct',df_grp1_demo)
print_mean_sd('Visual_correct',df_grp1_demo)
print_mean_sd('LNS',df_grp1_demo)
print_mean_sd('Spatial Span',df_grp1_demo)

print('SCZ grp       mean        std')
print_mean_sd('IQ',df_grp2_demo)
print_mean_sd('TMT_A',df_grp2_demo)
print_mean_sd('TMT_B',df_grp2_demo)
print_mean_sd('Verbal_correct',df_grp2_demo)
print_mean_sd('Visual_correct',df_grp2_demo)
print_mean_sd('LNS',df_grp2_demo)
print_mean_sd('Spatial Span',df_grp2_demo)

print('SCZ grp       mean        std')
print_mean_sd('B_PANSS_T',df_grp2_demo)
print_mean_sd('B_PANSS_P',df_grp2_demo)
print_mean_sd('B_PANSS_N',df_grp2_demo)
print_mean_sd('B_PANSS_G',df_grp2_demo)



### figures


sns.regplot(G1_mean_data1, df_grp1_TMT_A)
sns.regplot(G1_mean_data2, df_grp2_TMT_A)
plt.xlabel('1st gradient (10$^{-5}$)')
plt.ylabel('TMT-A time')
plt.xlim([-0.00003, 0.000055])
plt.xticks(ticks=[-0.00002, 0, 0.00002, 0.00004])

sns.regplot(G1_mean_data1, df_grp1_TMT_B)
sns.regplot(G1_mean_data2, df_grp2_TMT_B)
plt.xlabel('1st gradient (10$^{-5}$)')
plt.ylabel('TMT-B time')
plt.xlim([-0.00003, 0.00005])
plt.ylim([0, 400])
# plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
plt.xticks(ticks=[-0.00002, 0, 0.00002, 0.00004])

plt.savefig(f'{out_dir}/G1_TMT_B_Correlation.png') # needs fixing

