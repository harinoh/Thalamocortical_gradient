%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Construct stepwise connectivity matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sfc()
addpath(genpath('/Volume/CCNC/harin_oh/1_thalamocortical/code/SFC'));

disp(['## stepwise functinal connectivity - processing', newline]);
basepath = '/Volume/CCNC/harin_oh/1_thalamocortical/';
outpath = '/Volume/CCNC/harin_oh/1_thalamocortical/SFC_result/';
sublist = importdata([basepath, 'SCZOHC/Subject_list_SCZOHC.txt']);
Nsub = length(sublist);
Nstep = 200;

seed_name = 'whlThal'; % ### change the seed name
% seed = strsplit(seed_name,'_');
% seed = seed{1};
thal_idx = load([outpath,seed_name,'_seed_idx_from_ctxthal.mat']).indx_thal_ctxthal; % #### change the seed, 'ctxthal_',
if ~exist([outpath,'2_sfc/',seed_name],'dir') % creating folder
    mkdir([outpath,'2_sfc/',seed_name])
end

%% 2) Construct SFC matrix
disp(['## stepwise functinal connectivity - ',seed_name, newline]);
%sfc_all = cell(Nsub,1);
for sidx = 1 : Nsub
    subID = sublist{sidx};
    sfc_filename = [outpath, '2_sfc/',seed_name '/', subID,'_sfc.mat'];     % require around 69GB CPU memory 
    if ~exist(sfc_filename, 'file')
        disp(['subject = ', num2str(subID)])
        conn = load([outpath, '0_Connectivity_matrix_thal2thalctx/',subID, '_WholeCTX_coef_FDR_corr.mat']).whl_thal_conn_matrix_correct;
        %binconn = load([outpath, '1_binconn/', subID, '_binconn.mat']).mask; % #### Change the end depending the seed
        %seed_conn = conn(seed_idx,:);
        sfc = findwalks_mod(conn, Nstep, thal_idx);
%         binconn = load([outpath, '1_binconn/sym_ver/', subID, '_binconn.mat']).binconn;
%         sfc = compute_sfc(binconn, Nstep);
        
        save(sfc_filename, 'sfc', '-v7.3');
        clear sfc;
        disp(['Subject ', subID, ' sfc successfully saved']);
    else
        disp(['Subject ', subID, ' sfc file exist']);
    end
end
end
