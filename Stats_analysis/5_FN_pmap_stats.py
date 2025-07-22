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

ROI= 'PFCSS_Mot_cerebell68' # PFCSS_Mot  PFCSS_Mot_cerebell68  WholeBrain_cerebell68
ref_group='SCZOHC'
par_num=1000
opt=f'_FWHM4_Parcel{par_num}' 
sub_group1 = 'OHC' # NOR
sub_group2 = 'SCZ' # patient
home_dir = '/Volumes/CCNC/harin_oh/1_Thalamocortical'
pmap_path = f'{home_dir}/1_Gradient_result/5_pmap_alignment/{ROI}{opt}'
file_name_group1 = f'{home_dir}/0_raw_data/Subject_list_{sub_group1}.txt'
file_name_group2 = f'{home_dir}/0_raw_data/Subject_list_{sub_group2}.txt'
out_dir = f'{home_dir}/1_Gradient_result/6_2_FN/SCZOHC/{ROI}{opt}'

brainstat_path = os.path.abspath(os.path.join(f'{home_dir}/codes_Final/BrainStat-master'))

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

# Load Out roi atlas
roi_atlas = nib.load(f'{home_dir}/ROI_atlas/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm_{ROI}.nii.gz')
roi_data = roi_atlas.get_fdata()

nVoxels = np.prod(roi_data.shape)
roi_data = np.reshape(roi_data,nVoxels)
par_val = np.unique(roi_data[roi_data !=0]) 

# load cerebellar 
cerebell_atlas = nib.load(f'{home_dir}/ROI_atlas/Buckner_2011/atl-Buckner7_space-MNI_dseg_flirt_add1000.nii.gz').get_fdata()
cerebell_atlas = np.reshape(cerebell_atlas,nVoxels)

cerebell68_val = [val for val in par_val if val >= 1001]
cerebell_lab = []
for c in cerebell68_val:
    c_ind = np.where(roi_data == c)[0]
    fn_val = cerebell_atlas[c_ind]
    cerebell_fn_val, cerebell_fn_count = np.unique(fn_val[fn_val !=0], return_counts=True)
    max_ind = np.argmax(cerebell_fn_count)
    cer68_ind = cerebell_fn_val[max_ind]
    cerebell_lab.append(cer68_ind)

# read functional network label
fn_label = pd.read_csv(f'{home_dir}/ROI_atlas/Schaefer2018_fsleyes_lut/Schaefer2018_1000Parcels_7Networks_order.lut', sep='\s+', header=None, names=['Parcel', 'R', 'G', 'B', 'Network'])
fn_label['Network_Name'] = fn_label['Network'].str.extract(r'7Networks_[LR]H_([A-Za-z]+)_')
# Renaming the labels
label_mapping = {'Vis' : 'VIS', 'SomMot' : 'SMN', 'DorsAttn':'DAN', 'SalVentAttn':'VAN', 'Limbic':'LIM','Cont':'FPN', 'Default': 'DMN'}
fn_label['Network_Name'] = fn_label['Network_Name'].replace(label_mapping)
# functional network label corresponding to each parcel
network_parcel = {network: fn_label[fn_label['Network_Name'] == network]['Parcel'].tolist() for network in fn_label['Network_Name'].unique()}

# loading cerebellar fn labels
cereb68_label_mapping = {1001:'VIS',1002:'SMN',1003:'DAN',1004:'VAN',1005:'LIM',1006:'FPN',1007:'DMN'}
cerebell68_fn_label = [cereb68_label_mapping[val] if val in cereb68_label_mapping else val for val in cerebell_lab]

# designating fn label for each parcel (ROI)
cortex_par_val = [val for val in par_val if val <= 1000]
par_label = []
for p in cortex_par_val:
    fn = fn_label.loc[fn_label['Parcel'] == p, 'Network_Name'].values[0]
    par_label.append(fn)

par_label = np.concatenate((par_label,cerebell68_fn_label)) # if cerebellum is used within ROI
par_label_fn, counts = np.unique(par_label, return_counts=True)
desired_order = ['SMN', 'DAN', 'VAN', 'LIM', 'FPN', 'DMN']
reordered_fn_label = np.array([s for s in desired_order if s in par_label_fn]).tolist()



###############################
#### 1. Extract cortical region
def pmap_extract(G_data, roi_data, subject_list):
    G_data_roi = {}
    G_data_roi_list = []
    for sub_id, sub in enumerate(subject_list):
        sub_data = G_data[sub]
        sub_data = np.reshape(sub_data,nVoxels)
        mean_par = np.zeros(len(par_val))
        for ind, val in enumerate(par_val):
            roi_indx = np.where(roi_data == val)[0]
            sub_roi = sub_data[roi_indx]
            mean_par[ind] = np.mean(sub_roi)
        G_data_roi[sub] = mean_par
        G_data_roi_list.append(mean_par)
            
    return G_data_roi, G_data_roi_list

# L-thalamus
G1_L_data1_pmap, G1_L_data1_pmap_list = pmap_extract(G1_L_data1, roi_data, subject_list_group1)
G1_L_data2_pmap, G1_L_data2_pmap_list = pmap_extract(G1_L_data2, roi_data, subject_list_group2)

# R-thalamus
G1_R_data1_pmap, G1_R_data1_pmap_list = pmap_extract(G1_R_data1, roi_data, subject_list_group1)
G1_R_data2_pmap, G1_R_data2_pmap_list = pmap_extract(G1_R_data2, roi_data, subject_list_group2)


#### 1-2.Normalization
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

def regress_glm(G_data1_list, G_data2_list):
    G1_data1_arr = np.array(G_data1_list)
    G1_data2_arr = np.array(G_data2_list)
    all_array = np.row_stack((G1_data1_arr, G1_data2_arr))
    
    # Covariate + 1 designated to X
    X = np.column_stack((np.ones_like(demo['age']), demo['Group'], demo['age'])) # demo['Group'], , demo['sex_num']
    model = LinearRegression()
    model.fit(X, all_array)
    
    resid = all_array - model.predict(X)
    data1_arr = resid[:len(G_data1_list),:]
    data2_arr = resid[len(G_data1_list):,:]
    
    data1_dic = {sub: data1_arr[i] for i, sub in enumerate(subject_list_group1)}
    data2_dic = {sub: data2_arr[i] for i, sub in enumerate(subject_list_group2)}
    
    return resid, data1_dic, data2_dic

G1_L_resid, G1_L_data1_resid, G1_L_data2_resid = regress_glm(G1_L_data1_pmap_list, G1_L_data2_pmap_list)
G1_R_resid, G1_R_data1_resid, G1_R_data2_resid = regress_glm(G1_R_data1_pmap_list, G1_R_data2_pmap_list)



###############################
#### 3. Extract gradient values by each functional network
def fn_extract(G_data, subject_list):
    fn_sizes = {label: counts[np.where(par_label_fn == label)[0][0]] for label in par_label_fn}
    fn_dict =  {key: np.zeros((len(subject_list),fn_sizes[key])) for key in par_label_fn}
    fn_dict_avg = {key: np.zeros((len(subject_list))) for key in par_label_fn}
    for sub_id, sub in enumerate(subject_list):
        sub_data = G_data[sub]
        fn_counters = {key: 0 for key in par_label_fn}
        for ind in range(0, len(par_val)):
            sub_par_val = sub_data[ind]
            fn_name = par_label[ind]
            # append to final variable
            fn_dict[fn_name][sub_id, fn_counters[fn_name]] = sub_par_val
            fn_counters[fn_name] += 1
        for name in par_label_fn:
            sub_fn = fn_dict[name][sub_id]
            fn_dict_avg[name][sub_id] = np.mean(sub_fn, axis=0)
    return fn_dict, fn_dict_avg

## w/ regression
# group1
G1_L_data1_fn, G1_L_data1_fn_avg = fn_extract(G1_L_data1_resid, subject_list_group1)
G1_R_data1_fn, G1_R_data1_fn_avg = fn_extract(G1_R_data1_resid, subject_list_group1)

# group2
G1_L_data2_fn, G1_L_data2_fn_avg = fn_extract(G1_L_data2_resid, subject_list_group2)
G1_R_data2_fn, G1_R_data2_fn_avg = fn_extract(G1_R_data2_resid, subject_list_group2)

## w/o regression
# group1
G1_L_data1_fn, G1_L_data1_fn_avg = fn_extract(G1_L_data1_pmap, subject_list_group1)
G1_R_data1_fn, G1_R_data1_fn_avg = fn_extract(G1_R_data1_pmap, subject_list_group1)

# group2
G1_L_data2_fn, G1_L_data2_fn_avg = fn_extract(G1_L_data2_pmap, subject_list_group2)
G1_R_data2_fn, G1_R_data2_fn_avg = fn_extract(G1_R_data2_pmap, subject_list_group2)

del G1_L_data1, G1_R_data1, G1_L_data2, G1_R_data2,



# reorder and format into list variable

def reorder_dict_2list(G_data_avg, desired_order):
    reorder_dict = {key: G_data_avg[key] for key in desired_order if key in G_data_avg}
    val_list = list(reorder_dict.values())
    return val_list

# fn version
# group1
G1_L_data1_fn_list = reorder_dict_2list(G1_L_data1_fn, desired_order)
G1_R_data1_fn_list = reorder_dict_2list(G1_R_data1_fn, desired_order)

# group2
G1_L_data2_fn_list = reorder_dict_2list(G1_L_data2_fn, desired_order)
G1_R_data2_fn_list = reorder_dict_2list(G1_R_data2_fn, desired_order)

# fn average version
# group1
G1_L_data1_fn_avg_list = reorder_dict_2list(G1_L_data1_fn_avg, desired_order)
G1_R_data1_fn_avg_list = reorder_dict_2list(G1_R_data1_fn_avg, desired_order)

# group2
G1_L_data2_fn_avg_list = reorder_dict_2list(G1_L_data2_fn_avg, desired_order)
G1_R_data2_fn_avg_list = reorder_dict_2list(G1_R_data2_fn_avg, desired_order)




###############################
### 3. Statistical testing
from scipy import stats
from statsmodels.stats.multitest import multipletests

## t-test version
def ttest_FDR_correction(G_data1, G_data2, g):
    print(f'For {g} gradient')
    print(f'Statistical result when {ROI} ROI used')
    net_pVal_G = []
    net_tstat_G = []
    for n in range(0, len(par_label_fn)):
        x = G_data1[n]
        y = G_data2[n]
        
        t_stat, p_val = stats.ttest_ind(x, y, equal_var=False)
        net_pVal_G.append(p_val)
        net_tstat_G.append(t_stat)
        
        print(f'For network {n} in {g} gradient :')
        print('t-stat     p-value', )
        print('{:.4f}'.format(t_stat), '    ', '{:.5f}'.format(p_val))
    
    # Perform Bonferroni corretion
    print('Bonferroni-corrected p-value:')
    Bonf_corr_p_val = multipletests(net_pVal_G, method='bonferroni', alpha = 0.05)[1]
    for i, fn in enumerate(par_label_fn):
        print(f'Network {fn}: ', '{:.5f}'.format(Bonf_corr_p_val[i]))
    
    FDR_corrected_p_val = multipletests(net_pVal_G, method='fdr_bh', alpha = 0.05)[1]
    print('FDR-corrected (fdr_bh) p-value:')
    for i, fn in enumerate(par_label_fn):
        print(f'Network {fn}: ', '{:.5f}'.format(FDR_corrected_p_val[i]))    
    return net_pVal_G, FDR_corrected_p_val #Bonf_corr_p_val,

net_L_pVal_G1, G1_L_FDR_p_val = ttest_FDR_correction(G1_L_data1_fn_avg_list, G1_L_data2_fn_avg_list, '1st')


## Surfstat version
from brainstat.stats.SLM import SLM
from brainstat.stats import terms

def brainstat_avg_ttest(G_data1_fn_list, G_data2_fn_list):
    term_intercept = terms.FixedEffect(1, names='intercept')
    term_age = terms.FixedEffect(demo['age'], names='age')
    term_group = terms.FixedEffect(demo['Group'], names='group')
    term_sex = terms.FixedEffect(demo['sex_num'], names='sex')
    #term_fdrms = terms.FixedEffect(concat_mean_fdrms)
    
    group = demo['Group']
    model = term_intercept + term_group + term_age #+ term_sex # + term_fdrms 
    slm_age = SLM(model, -group, correction='fdr')
    
    G_data1_fn = np.array(G_data1_fn_list).T
    G_data2_fn = np.array(G_data2_fn_list).T
    G_concat = np.concatenate((G_data1_fn, G_data2_fn), axis=0)

    slm_age.fit(G_concat)
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
    
    for n in range(len(par_label_fn)):
        print( desired_order[n],':', '{:.5f}'.format(ttest[:,n][0]),'    ','{:.5f}'.format(pval[n,:][0]),'     ', '{:.5f}'.format(qval[n]))
    
    return ttest, qval, fdr_pval_idx, pval

# Avg version
G1_L_ttest, G1_L_qval, G1_fdr_pval_idx, G1_L_pval = brainstat_avg_ttest(G1_L_data1_fn_avg_list, G1_L_data2_fn_avg_list)
G1_R_ttest, G1_R_qval, G1_fdr_pval_idx, G1_R_pval = brainstat_avg_ttest(G1_R_data1_fn_avg_list, G1_R_data2_fn_avg_list)


###############################
##### Figure
def avg_dict(G_data_dict):
    fn_sizes = {label: counts[np.where(par_label_fn == label)[0][0]] for label in par_label_fn}
    avg_fn_dict = {key: np.zeros((fn_sizes[key])) for key in par_label_fn}
    for fn in par_label_fn:
        G_data_dict_fn =  G_data_dict[fn]
        avg_fn = np.mean(G_data_dict_fn, axis=0)
        avg_fn_dict[fn] = avg_fn
    
    return avg_fn_dict

# L-thalamus
G1_L_data1_avg = avg_dict(G1_L_data1_fn)
G1_L_data2_avg = avg_dict(G1_L_data2_fn)

# R-thalamus
G1_R_data1_avg = avg_dict(G1_R_data1_fn)
G1_R_data2_avg = avg_dict(G1_R_data2_fn)



import matplotlib.pyplot as plt
import seaborn as sns

desired_order_plt = ['DAN', 'LIM', 'SMN', 'DMN', 'VAN', 'FPN']
#desired_order_plt = ['LIM', 'DMN', 'FPN', 'VAN', 'SMN', 'DAN']
#desired_order_plt = ['SMN', 'LIM', 'DAN', 'DMN', 'VAN', 'FPN']
reordered_fn_label = np.array([s for s in desired_order_plt if s in par_label_fn]).tolist()


def plt_densMap_combine(G1_L_data1_avg, G1_L_data2_avg, G, figsize=(8,0.8), share_x=True, xlim=(-100, 100), save_path=None):
    num_net = len(reordered_fn_label)
    fig, axs = plt.subplots(num_net, 1, figsize=(figsize[0], num_net*figsize[1]), sharex=share_x)
    colours = ['#138144','#d2f688','#4388B1','#CB4A60','#d63fff', '#E59036'] # DAN > LIM > SMN > DMN > VAN > FPN
    grp2_col = ['#a6f2c8','#fafde7','#b5d2e3', '#e8b0b9', '#e999ff',  '#f3cda5']
    line_col_grp2 = ['#00833D','#94d510','#0070B1', '#CD0022', '#c800ff',  '#E67700']
    # colours = ['#d2f688','#CB4A60', '#E59036','#d63fff','#4388B1', '#138144'] # LIM > DMN > FPN > VAN > SMN > DAN
    # grp2_col = ['#fafde7','#e8b0b9','#f3cda5','#e999ff','#b5d2e3', '#a6f2c8']
    # line_col_grp2 = ['#94d510', '#CD0022','#E67700', '#c800ff', '#0070B1', '#00833D']
    # colours = ['#4388B1','#d2f688','#138144','#CB4A60','#d63fff','#E59036'] #SMN > LIM > DAN > DMN > VAN > FPN
    # grp2_col = ['#b5d2e3', '#fafde7', '#a6f2c8', '#e8b0b9', '#e999ff', '#f3cda5']
    # line_col_grp2 = ['#0070B1', '#94d510', '#00833D', '#CD0022', '#c800ff', '#E67700']
    # colours = ['#d2f688','#CB4A60','#d63fff', '#E59036','#4388B1', '#138144'] # LIM > DMN > VAN > FPN > SMN > DAN
    # grp2_col = ['#fafde7','#e8b0b9','#e999ff','#f3cda5','#b5d2e3', '#a6f2c8']
    # line_col_grp2 = ['#94d510', '#CD0022','#c800ff','#E67700', '#0070B1', '#00833D']
    for i, net in enumerate(reordered_fn_label):
        sns.kdeplot(data=G1_L_data2_avg[net], ax=axs[i], color=grp2_col[i], shade=True, alpha=0.8)
        sns.kdeplot(data=G1_L_data1_avg[net], ax=axs[i], color=colours[i], shade=True, alpha=0.8)
        sns.kdeplot(data=G1_L_data2_avg[net], ax=axs[i], color=line_col_grp2[i], alpha=1)
      
        # removing the box around each figure
        axs[reordered_fn_label.index(net)].spines['top'].set_color('none')
        axs[reordered_fn_label.index(net)].spines['right'].set_color('none')
        axs[reordered_fn_label.index(net)].spines['left'].set_color('none')
        # only labeling x-axis at the bottom
        if i == 5:
            axs[reordered_fn_label.index(net)].set_xlabel(f'{G} density', labelpad=5)
        # rotating the y-axis label & aligning to the center
        axs[i].yaxis.label.set_rotation(0)
        axs[i].yaxis.label.set_horizontalalignment('center')
        axs[i].set_ylabel(f'{reordered_fn_label[i]}')
        # remove y_axis label and sticks
        axs[i].set_yticks([])
        axs[i].set_yticklabels([])
        # setting min max for x-axis
        axs[i].set_xlim(xlim[0], xlim[1])
        axs[i].set_xticks(xlim)
    plt.tight_layout()
    plt.show()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        

ext = 'Vox'  # Vox_norm
plt_densMap_combine(G1_L_data1_avg, G1_L_data2_avg, 'G1: L-thal projection', figsize=(8,0.8), xlim=(-1800, 2200), save_path=f'{out_dir}/Figure/FN_distribution_G1_L_thal_avg{ext}.png') # -1800, 2200
plt_densMap_combine(G1_R_data1_avg, G1_R_data2_avg, 'G1: R-thal projection', figsize=(8,0.8), xlim=(-1800, 2500), save_path=f'{out_dir}/Figure/FN_distribution_G1_R_thal_avg{ext}.png') # -1800, 2500



## violin plot
from matplotlib.lines import Line2D

def plt_violin_combine(G1_L_data1_avg, G1_L_data2_avg, G, save_path=None):
    hc_colors = ['#138144','#d2f688','#4388B1','#CB4A60','#d63fff', '#E59036']  # DAN > LIM > SMN > DMN > VAN > FPN
    scz_colors = ['#a6f2c8','#fafde7','#b5d2e3', '#e8b0b9', '#e999ff',  '#f3cda5']
    # hc_colors = ['#d2f688','#CB4A60', '#E59036','#d63fff','#4388B1', '#138144'] # LIM > DMN > FPN > VAN > SMN > DAN
    # scz_colors = ['#fafde7','#e8b0b9','#f3cda5','#e999ff','#b5d2e3', '#a6f2c8']

    sns.set(style="ticks")
    fig, ax = plt.subplots(figsize=(9, 0.8 * len(reordered_fn_label)))

    # Plot each network individually
    for i, net in enumerate(reordered_fn_label):
        data1 = G1_L_data1_avg[net]
        data2 = G1_L_data2_avg[net]
        df = pd.DataFrame({
            'Value': list(data1) + list(data2),
            'Group': ['HC'] * len(data1) + ['SCZ'] * len(data2),
            'Network': [net] * (len(data1) + len(data2))
        })

        sns.violinplot(
            data=df,
            x='Value',
            y='Network',
            hue='Group',
            split=True,
            inner='quart',
            gap=0.1,
            palette={'HC': hc_colors[i], 'SCZ': scz_colors[i]},
            ax=ax
        )
    
    for art in ax.collections:  # adjust number if needed
        art.set_edgecolor('#575757')

    ax.set(xlabel=f'{G} value', ylabel='Functional network')
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(axis='y', which='both', length=0)
    ax.set_xlim(-1500, 2300)
    ax.set_xticks((-1500, 2300))

    # === Custom legend with arrow markers ===
    arrow_up = Line2D([0], [0], marker=r'$\uparrow$', color='black', linestyle='None', markersize=8)
    baseline = Line2D([0], [0], color='black', lw=1)
    arrow_down = Line2D([0], [0], marker=r'$\downarrow$', color='black', linestyle='None', markersize=8)

    ax.legend(
        handles=[arrow_up, baseline, arrow_down],
        labels=['HC', '', 'SCZ'],
        title='Group',
        loc='upper left',
        bbox_to_anchor=(1.02, 1),
        borderaxespad=0.,
        frameon=True,
        handlelength=1.5,
        handletextpad=0.8,
        labelspacing=0.3
    )

    plt.tight_layout()
    plt.show()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
    
        
name='avg'
plt_violin_combine(G1_L_data1_avg, G1_L_data2_avg, 'G1: L-thal projection', save_path=f'{out_dir}/Figure/FN_Pmap_violinplot_G1_L_thal_{name}.png')
plt_violin_combine(G2_L_data1_avg, G2_L_data2_avg, 'G2: L-thal projection', save_path=f'{out_dir}/Figure/FN_Pmap_violinplot_G2_L_thal_{name}.png')

