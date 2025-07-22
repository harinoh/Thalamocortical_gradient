#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 21:50:05 2023

@author: harin_oh
"""


import nibabel as nib
import numpy as np
import pandas as pd
from sklearn.preprocessing import normalize
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os
import sys

group = 'SCZOHC'
roi= 'WholeCTX'
ref_group='SCZOHC'
opt='_FWHM4_IndivMul'
sub_group1 = 'OHC' # NOR
sub_group2 = 'SCZ' # patient
home_dir = '/Volume/CCNC/harin_oh/1_thalamocortical'
gradient_path = f'{home_dir}/Gradient_result/2_Alignment/{group}/{roi}_{ref_group}{opt}'
file_name_group1 = f'{home_dir}/{group}/Subject_list_{sub_group1}.txt'
file_name_group2 = f'{home_dir}/{group}/Subject_list_{sub_group2}.txt'
out_dir = f'{home_dir}/Gradient_result/3_2_FunctionalNetwork/{group}/{roi}_{ref_group}{opt}/Figure'
brainstat_path = os.path.abspath(os.path.join(f'{home_dir}/code/BrainStat-master'))

if brainstat_path not in sys.path:
    sys.path.append(brainstat_path)
    
# def function
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
        L_thal = nib.load(f'{gradient_path}/{sub}_Gordon_L_thal_222_aligned.cmaps.nii.gz')
        R_thal = nib.load(f'{gradient_path}/{sub}_Gordon_R_thal_222_aligned.cmaps.nii.gz')
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




##### 0. Preparation ###
# make directory if not available
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    print(f'Folder {out_dir} created.')
else:
    print(f'Folder {out_dir} already exists')


print('Loading fMRI images & thalamic atlas')
# read subject file
subject_list_group1 = read_subject_list(file_name_group1)
subject_list_group2 = read_subject_list(file_name_group2)
print(f'NOR group: {subject_list_group1}')
print(f'patient group: {subject_list_group2}')

# Loading thalamic gradient
print(f'Loading cmaps from {gradient_path}')
G1_L_data1, G2_L_data1, G3_L_data1, G1_R_data1, G2_R_data1, G3_R_data1 = nii_read(subject_list_group1)      # group NOR
G1_L_data2, G2_L_data2, G3_L_data2, G1_R_data2, G2_R_data2, G3_R_data2 = nii_read(subject_list_group2)      # group SCZ

# Loading demographic data
cols = ['Folder_name', 'Group','age', 'sex_num', 'IQ'] #, 'sex_num'
demo_grp1 = pd.read_excel(f'{home_dir}/BCS_ECT_OHC_rsfMRI_Demo.xlsx', sheet_name='baseline_n129', usecols=cols, nrows=69, index_col='Folder_name')
demo_grp2 = pd.read_excel(f'{home_dir}/BCS_ECT_OHC_rsfMRI_Demo.xlsx', sheet_name='baseline_n129', usecols=cols, skiprows=range(1,70), nrows=60, index_col='Folder_name')
demo = pd.concat([demo_grp1,demo_grp2])

# Load thalamic functional atlas
L_thal_atlas=nib.load(f'/Volume/CCNC_T_4/harin_oh/1_thalamocortical/ROI_atlas/Gordon_Parcels/Gordon_L_thal_222.nii.gz')
R_thal_atlas=nib.load(f'/Volume/CCNC_T_4/harin_oh/1_thalamocortical/ROI_atlas/Gordon_Parcels/Gordon_R_thal_222.nii.gz')
L_thal = L_thal_atlas.get_fdata()
R_thal = R_thal_atlas.get_fdata()
L_thal_indx = np.where(L_thal !=0)
R_thal_indx = np.where(R_thal !=0)

opt2 = '_PFCSS_Mot' #'_Yeo7_Network' #'' #'_PFCSS_Mot'
L_thal_atlas=f'/Volume/CCNC_T_4/harin_oh/1_thalamocortical/Gradient_result/3_2_FunctionalNetwork/Thalamic_FN_Gordon_L_thal_222_{group}_FWHM4{opt2}.nii.gz'
R_thal_atlas=f'/Volume/CCNC_T_4/harin_oh/1_thalamocortical/Gradient_result/3_2_FunctionalNetwork/Thalamic_FN_Gordon_R_thal_222_{group}_FWHM4{opt2}.nii.gz'
print(f'loading thalamus functional atlas: {L_thal_atlas}, {R_thal_atlas}')
L_thal_FN = nib.load(L_thal_atlas).get_fdata()
R_thal_FN = nib.load(R_thal_atlas).get_fdata()
L_thal_FN_extract = L_thal_FN[L_thal_indx]
R_thal_FN_extract = R_thal_FN[R_thal_indx]
L_thal_FN_indx = np.where(L_thal_FN_extract !=0)
R_thal_FN_indx = np.where(R_thal_FN_extract !=0)

# Load FDRMS file
def read_FDRMS(subject_list):
    fdrms = {}
    for i in subject_list:
        fdrms_dir = f'{home_dir}/{group}/{i}/preproc/{i}_MoCorr_FDRMS_val_MO'
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



###############################
### 1. extract BOLD signal value according to functional network atlas

### 1-1. extracting FN values
FN_val = np.unique(L_thal_FN)
FN_val_list = [int(value) for value in FN_val if value > 0]
# FN_val_list = [int(value) for value in FN_val if value > 0 and value != 5]



### 1-2. extract thalamic BOLD signal
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


del G1_L_data1, G1_R_data1, G2_L_data1, G2_R_data1, G3_L_data1, G3_R_data1, G1_L_data2, G1_R_data2, G2_L_data2, G2_R_data2, G3_L_data2, G3_R_data2


# if IQ is considered as covariate and IQ has nan values
# from sklearn.impute import SimpleImputer
# imputer = SimpleImputer(strategy='mean')  # You can use other strategies like 'median', 'most_frequent', etc.
# imputed_demo = pd.DataFrame(imputer.fit_transform(demo), columns=demo.columns)

### 1-3. Regress
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


### 1-4. extract values into FN
def FN_extract(G_L_data, G_R_data, L_thal_FN, R_thal_FN, subject_list):
    if isinstance(G_L_data, dict):
        G_L_data1_list = list(G_L_data.values())
        G_R_data1_list = list(G_R_data.values())
    elif isinstance(G_L_data, np.ndarray):
        G_L_data1_list = G_L_data
        G_R_data1_list = G_R_data
    G_data_FN = [[] for _ in range(len(FN_val_list))]
    for fn in FN_val_list:
        L_thal_indx_FN = np.where(L_thal_FN == fn)
        R_thal_indx_FN = np.where(R_thal_FN == fn)
        G_data_sub_thal = []
        for sub in range(0, len(subject_list)):
            sub_L_thal = G_L_data1_list[sub][L_thal_indx_FN]
            sub_R_thal = G_R_data1_list[sub][R_thal_indx_FN]
            sub_thal = np.concatenate((sub_L_thal,sub_R_thal))
            G_data_sub_thal.append(sub_thal)
        G_data_sub_thal = np.array(G_data_sub_thal)
        if fn >=6:
            G_data_FN[fn-3] = G_data_sub_thal
        else:
            G_data_FN[fn-2] = G_data_sub_thal
        
    return G_data_FN

# w/ regression
G1_data1_FN = FN_extract(G1_L_resid_data1, G1_R_resid_data1, L_thal_FN_extract, R_thal_FN_extract, subject_list_group1)
G1_data2_FN = FN_extract(G1_L_resid_data2, G1_R_resid_data2, L_thal_FN_extract, R_thal_FN_extract, subject_list_group2)

# w/o regression
G1_data1_FN = FN_extract(G1_L_data1_thal, G1_R_data1_thal, L_thal_FN_extract, R_thal_FN_extract, subject_list_group1)
G1_data2_FN = FN_extract(G1_L_data2_thal, G1_R_data2_thal, L_thal_FN_extract, R_thal_FN_extract, subject_list_group2)


 
###############################
#### 2. Average voxel by each functional network
def mean_2_array(gradient, FN_val_list):
    net_avg_grp_G1 = []
    for n in range(0, len(FN_val_list)):
        net1_G = gradient[n]
        avg = np.mean(net1_G,axis=1)
        net_avg_grp_G1.append(avg)
    net_avg_grp_G1 = np.array(net_avg_grp_G1)
    return net_avg_grp_G1

net_avg_G1_data1 = mean_2_array(G1_data1_FN, FN_val_list)
net_avg_G1_data2 = mean_2_array(G1_data2_FN, FN_val_list)



###############################
### 3. Statistcal significance between two groups

# t-test
def ttest_FDR_correction(FN_val_list, array_G_data1, array_G_data2, roi, g):
    print(f'For {g} gradient')
    print(f'Statistical result when {roi} ROI used')
    net_pVal_G = []
    net_tstat_G = []
    for n in range(0, len(FN_val_list)):
        x = array_G_data1[n]
        y = array_G_data2[n]
        
        t_stat, p_val = stats.ttest_ind(x, y, equal_var=False)
        net_pVal_G.append(p_val)
        net_tstat_G.append(t_stat)
        
        print(f'For network{n+1} in {g} gradient :')
        print('t-stat     p-value', )
        print('{:.4f}'.format(t_stat), '    ', '{:.5f}'.format(p_val))
    
    # Perform Bonferroni corretion
    print('Bonferroni-corrected p-value:')
    Bonf_corr_p_val = multipletests(net_pVal_G, method='bonferroni', alpha = 0.05)[1]
    for i in range(0, len(FN_val_list)):
        print(f'Network {i}: ', '{:.5f}'.format(Bonf_corr_p_val[i]))
    
    FDR_corrected_p_val = multipletests(net_pVal_G, method='fdr_bh', alpha = 0.05)[1]
    print('FDR-corrected (fdr_bh) p-value:')
    for i in range(0, len(FN_val_list)):
        print(f'Network {i}: ', '{:.5f}'.format(FDR_corrected_p_val[i]))
        
    return net_pVal_G, Bonf_corr_p_val, FDR_corrected_p_val #Bonf_corr_p_val,


# average ver
net_pVal_G1, G1_Bonferroni_p_val, G1_FDR_p_val = ttest_FDR_correction(FN_val_list, net_avg_G1_data1, net_avg_G1_data2, roi, '1st') 
net_pVal_G2, G2_Bonferroni_p_val, G2_FDR_p_val = ttest_FDR_correction(FN_val_list, net_avg_G2_data1, net_avg_G2_data2, roi, '2nd')  
net_pVal_G3, G3_Bonferroni_p_val, G3_FDR_p_val = ttest_FDR_correction(FN_val_list, net_avg_G3_data1, net_avg_G3_data2, roi, '3rd') 


# ANCOVA
import statsmodels.formula.api as smf

def ANCOVA_construct(array_G_data1, array_G_data2, roi, demo, G):
    print(f'Statistical result when {roi} ROI used')
    f_stat = []
    p_val = []
    for n in range(0,len(FN_val_list)):
        grp1 = array_G_data1[n]
        grp2 = array_G_data2[n]

        df_grp1 = pd.DataFrame(grp1,
                            index=subject_list_group1,
                            columns=['outcome'])
        df_grp2 = pd.DataFrame(grp2, 
                            index=subject_list_group2, 
                            columns=['outcome'])
        df_all = pd.concat([df_grp1, df_grp2])
        data = pd.concat([df_all, demo], axis=1)
        
        # ANCOVA model fitting
        model = smf.ols('outcome ~ Group + age', data=data).fit() # + meanFDRMS + sex_num
        f_stat.append(model.fvalue)
        p_val.append(model.f_pvalue)
    
    print('For', G)
    Bonf_corr_p_val = multipletests(p_val, method='bonferroni', alpha = 0.05)[1]
    FDR_corr_p_val = multipletests(p_val, method='fdr_bh', alpha = 0.05)[1]
    for x in range(0,len(FN_val_list)):
        print('uncorrected p-val:', p_val[x])
        print('Bonf_corr voxels:' , Bonf_corr_p_val[x])
        print('FDR_corr voxels:', FDR_corr_p_val[x])
        print('                ')
    
    return f_stat, p_val

f_stat_G1, p_val_G1 = ANCOVA_construct(net_avg_G1_data1, net_avg_G1_data2, roi, demo, 'G1')
f_stat_G2, p_val_G2 = ANCOVA_construct(net_avg_G2_data1, net_avg_G2_data2, roi, demo, 'G2')
f_stat_G3, p_val_G3 = ANCOVA_construct(net_avg_G3_data1, net_avg_G3_data2, roi, demo, 'G3')



# Surfstat version
from brainstat.stats.SLM import SLM
from brainstat.stats import terms


def brainstat_ttest(demo, data1_fn_list, data2_fn_list, roi):
    print(f'Statistical result when {roi} ROI used')
    
    fn_ttest = []
    fn_qval = []
    # fn_pval = []
    for n in range(0,len(FN_val_list)):
        term_intercept = terms.FixedEffect(1, names='intercept')
        term_age = terms.FixedEffect(demo['age'], names='age')
        term_group = terms.FixedEffect(demo['Group'], names='group')
        #term_sex = terms.FixedEffect(demo['sex_num'], names='sex')
        #term_fdrms = terms.FixedEffect(concat_mean_fdrms)
    
        group = demo['Group']
        model = term_intercept + term_group + term_age #+ term_sex + term_fdrms 
        slm_age = SLM(model, -group, correction='fdr')
        
        # fn_data1 = data1_fn_list[n]
        # fn_data2 = data2_fn_list[n]
        fn_data1 = np.mean(data1_fn_list[n], axis=1)
        fn_data2 = np.mean(data2_fn_list[n], axis=1)
        fn_concat = np.concatenate((fn_data1,fn_data2), axis=0)
        fn_concat_2d = fn_concat.reshape(-1,1)
        # fn_concat = np.array(fn_concat)

        slm_age.fit(fn_concat_2d)
        ttest = slm_age.t[0]
        qval = slm_age.Q[0]
        
        fn_ttest.append(ttest)
        fn_qval.append(qval)
        
        # if slm_age.P is None:
        #     print("P-values are not calculated by default. Attempting manual calculation...")

        #     # Assuming T-statistics are available in the fitted model
        #     if hasattr(slm_age, 't'):
        #         from scipy.stats import t

        #         df = slm_age.df  # Degrees of freedom
        #         t_values = slm_age.t   # T-statistics

        #         # Calculate p-values
        #         p_values = 2 * (1 - t.cdf(np.abs(t_values), df))
        #         slm_age.P = p_values
        #         fn_pval.append(np.transpose(p_values)[0])
        #     else:
        #         print("T-statistics not available. Cannot calculate P-values manually.")
        # else:
        #     print("P-values are already calculated.")
    
    FDR_corr_p_val = multipletests(fn_qval, method='fdr_bh', alpha = 0.05)[1]
    FN_label = ['SMN', 'DAN', 'VAN', 'FPN', 'DMN', 'VIS']
    for n in range(0, len(FN_val_list)):
        print('For functional networks', FN_label[n])
        print('t stat:', '{:.5f}'.format(float(fn_ttest[n])))
        print('p value:', '{:.5f}'.format(float(fn_qval[n])))
        print('FDR-corrected p-value:', '{:.5f}'.format(float(FDR_corr_p_val[n])))
    
    return fn_ttest, fn_qval, FDR_corr_p_val

G1_ttest, G1_qval, G1_pval = brainstat_ttest(demo, G1_data1_FN, G1_data2_FN, roi)
G2_ttest, G2_qval, G2_pval = brainstat_ttest(demo, G2_data1_FN, G2_data2_FN, roi)




###############################
### 4. Construct density map for each thalamic functional network

# flat version
def flat_array(grp1_G_FN_array, grp2_G_FN_array):
    grp1_G_FN_flat = [[] for _ in range(len(FN_val_list))]
    grp2_G_FN_flat = [[] for _ in range(len(FN_val_list))]
    for fn in range(0, len(FN_val_list)):
        grp1_fn_arr = grp1_G_FN_array[fn]
        grp2_fn_arr = grp2_G_FN_array[fn]
        grp1_G_FN_flat[fn] = grp1_fn_arr.flatten()
        grp2_G_FN_flat[fn] = grp2_fn_arr.flatten()
    return grp1_G_FN_flat, grp2_G_FN_flat

grp1_G1_FN_flat, grp2_G1_FN_flat = flat_array(G1_data1_FN, G1_data2_FN)




def plt_densMap_combine(FN_val_list, net_gradient_grp1, net_gradient_grp2, G, figsize=(8,0.8), share_x=True, xlim=(-0.1, 0.1), save_path=None):
    num_net = len(FN_val_list)
    fig, axs = plt.subplots(num_net, 1, figsize=(figsize[0], num_net*figsize[1]), sharex=share_x)
    FN_val_list_name={0: 'VIS', 1: 'SMN', 2: 'DAN', 3: 'VAN', 4: 'FPN', 5: 'DMN'}
    colours = ['#703996', '#4388B1', '#138144', '#d63fff', '#E59036', '#CB4A60']
    grp2_col = ['#d0b5e3', '#b5d2e3', '#a6f2c8', '#e999ff', '#f3cda5', '#e8b0b9']
    line_col_grp2 = ['#5F0C97', '#0070B1', '#00833D', '#c800ff', '#E67700', '#CD0022']
    for i, net in enumerate(FN_val_list):
        FN_lab = FN_val_list_name[i]
        sns.kdeplot(data=net_gradient_grp2[i], ax=axs[i], color=grp2_col[i], shade=True, alpha=0.8)
        sns.kdeplot(data=net_gradient_grp1[i], ax=axs[i], color=colours[i], shade=True, alpha=0.6)
        sns.kdeplot(data=net_gradient_grp2[i], ax=axs[i], color=line_col_grp2[i], alpha=1)
      
        # removing the box around each figure
        axs[FN_val_list.index(net)].spines['top'].set_color('none')
        axs[FN_val_list.index(net)].spines['right'].set_color('none')
        axs[FN_val_list.index(net)].spines['left'].set_color('none')
        # only labeling x-axis at the bottom
        if net == 7:
            axs[FN_val_list.index(net)].set_xlabel(f'{G} Density', labelpad=5)
        # rotating the y-axis label & aligning to the center
        axs[i].yaxis.label.set_rotation(0)
        axs[i].yaxis.label.set_horizontalalignment('center')
        axs[i].set_ylabel(f'{FN_lab}')
        # remove y_axis label and sticks
        axs[i].set_yticks([])
        axs[i].set_yticklabels([])
        # setting min max for x-axis
        axs[i].set_xlim(xlim[0], xlim[1])
        axs[i].set_xticks(xlim)
    plt.tight_layout()
    plt.show()
    if save_path:
        fig.savefig(save_path, dpi=100, bbox_inches='tight')


norm_stat = 'flat_ageCo_fromYeo7FN_atlas'
plt_densMap_combine(FN_val_list, grp1_G1_FN_flat, grp2_G1_FN_flat, 'G1', figsize=(8,0.8), share_x=True, xlim=(-0.008,0.006), save_path=f'{out_dir}/{group}_{roi}_G1_densityMap_{norm_stat}.png')


# Density map of each group
def plt_densMap(FN_val_list, net_gradient, G, figsize=(8,4), share_x=True, xlim=(-0.1, 0.1), save_path=None):
    num_net = len(FN_val_list)
    fig, axs = plt.subplots(num_net, 1, figsize=(figsize[0], num_net*figsize[1]), sharex=share_x)
    # colours = ['#703996', '#4388B1', '#138144', '#d63fff', '#E59036', '#CB4A60']  # for 6 FN atlas
    # FN_val_list_name = {0: 'VIS', 1: 'SMN', 2: 'DAN', 3: 'VAN', 4: 'FPN', 5: 'DMN'} # for 6 FN atlas
    colours = ['#4388B1', '#138144', '#d63fff', '#E59036', '#CB4A60'] # for 5 FN atlas
    FN_val_list_name = {0: 'SMN', 1: 'DAN', 2: 'VAN', 3: 'FPN', 4: 'DMN'} # for 5 FN atlas
    for i, net in enumerate(FN_val_list):
        FN_lab = FN_val_list_name[i]
        sns.kdeplot(net_gradient[i], ax=axs[i], color=colours[i], shade=True)  # for SCZ grp
        # sns.kdeplot(net_gradient[i], ax=axs[i], color=colours[i], shade=True, alpha=0.8) # for NOR grp
        # sns.kdeplot(net_gradient[i], ax=axs[i], color=colours[i], alpha=1) # for NOR grp
        # # removing the box around each figure
        axs[FN_val_list.index(net)].spines['top'].set_color('none')
        axs[FN_val_list.index(net)].spines['right'].set_color('none')
        axs[FN_val_list.index(net)].spines['left'].set_color('none')
        # only labeling x-axis at the bottom
        if net == 7:
            axs[FN_val_list.index(net)].set_xlabel(f'{G} Density', labelpad=10)
        # rotating the y-axis label & aligning to the center
        axs[i].yaxis.label.set_rotation(0)
        axs[i].yaxis.label.set_horizontalalignment('center')
        axs[i].set_ylabel(f'{FN_lab}')
        axs[i].set_yticks([])
        axs[i].set_yticklabels([])
        axs[i].set_xlim(xlim[0], xlim[1])
        axs[i].set_xticks(xlim)
    plt.tight_layout()
    plt.show()
    if save_path:
        fig.savefig(save_path, dpi=100, bbox_inches='tight')
        
# group distribution
add ='_flat_ageCo'
plt_densMap(FN_val_list, grp1_G1_FN_flat, 'NOR G1', figsize=(8,0.8), share_x=True, xlim=(-0.008,0.0055), save_path=f'{out_dir}/{sub_group1}_{roi}_G1_densityMap{add}.png')
plt_densMap(FN_val_list, grp1_G2_FN_flat, 'NOR G2', figsize=(8,0.8), share_x=True, xlim=(-0.0105,0.0105), save_path=f'{out_dir}/{sub_group1}_{roi}_G2_densityMap{add}.png')

# Patient group
plt_densMap(FN_val_list, grp2_G1_FN_flat, 'SCZ G1', figsize=(8,0.8), share_x=True, xlim=(-0.008,0.0055), save_path=f'{out_dir}/{sub_group2}_{roi}_G1_densityMap{add}.png')
plt_densMap(FN_val_list, grp2_G2_FN_flat, 'SCZ G2', figsize=(8,0.8), share_x=True, xlim=(-0.0105,0.0105), save_path=f'{out_dir}/{sub_group2}_{roi}_G2_densityMap{add}.png')


## creating box plot
def plt_boxplot_diff(G_L_data1, G_L_data2, g, figsize=(5,5), ylim=(-0.013, 0.013)):
    data_concat = {'Value': np.concatenate((G_L_data1, G_L_data2)),
                    'Group': ['ROI0']*len(G_L_data1) + ['ROI1']*len(G_L_data2)}
    df = pd.DataFrame(data_concat)
    fig, axs = plt.subplots(figsize=(figsize[0], figsize[1]))
    sns.set(style="whitegrid")
    sns.boxplot(x='Group', y='Value', data=df, palette="Set2")
    sns.stripplot(x='Group', y='Value', data=df)
    plt.xlabel(g)
    axs.spines['top'].set_color('none')
    axs.spines['right'].set_color('none')
    axs.set_ylim(ylim[0], ylim[1])
    plt.show() 
