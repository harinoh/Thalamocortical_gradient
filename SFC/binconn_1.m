%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Construct stepwise connectivity matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stepwise_fc()
addpath(genpath('/Volume/CCNC/harin_oh/1_thalamocortical/code/SFC'));

basepath = '/Volume/CCNC/harin_oh/1_thalamocortical/';
outpath = '/Volume/CCNC/harin_oh/1_thalamocortical/SFC_result/';
sublist = importdata([basepath, 'SCZOHC/Subject_list_SCZOHC.txt']);
Nsub = length(sublist);
Nroi = 56447; % FO: 405, Int: 1357, HO: 630
Nstep = 10;

%% 1) Binarize connectivity matrix with threshold 5
disp(['## Binarization of connectivity - processing', newline]);

% load seed index
FO_seed = load([outpath, 'FO_seed_region.mat']).seed_idx_SMN;
Int_seed = load([outpath,'Int_seed_region.mat']).seed_idx_int;
HO_seed = load([outpath,'HO_seed_region.mat']).seed_idx_HO;

binconn_all = cell(Nsub,1);
parfor k = 1 : numel(sublist)
    sublist_k = sublist{k};
    sidx = sublist_k;
    binconn_filepath = [outpath, '1_binconn/', sidx];
   
    if ~exist([outpath,'1_binconn/', sidx, '_FO_binconn.mat'],'file')
        disp(['subject = ', num2str(sidx)])
        conn = load([outpath, '0_Connectivity_matrix_thal2thalctx/', sidx, '_WholeCTX_coef_FDR_corr.mat']).whl_thal_conn_matrix_correct;
        
        %extracting seed regions
        FO_conn = conn(FO_seed, :);
        Int_conn = conn(Int_seed, :);
        HO_conn = conn(HO_seed,:);
        
        % binarize
        FO_binconn = binarize_conn(FO_conn);
        Int_binconn = binarize_conn(Int_conn);
        HO_binconn = binarize_conn(HO_conn);
        binconn_all{k} = struct('subID', sidx, 'binconn_filepath', binconn_filepath, 'FO_binconn', FO_binconn, 'Int_binconn', Int_binconn, 'HO_binconn', HO_binconn);
        
        disp(['Subject ', sidx, ' binnary matrix construction completed']);
    else
        disp(['Subject ', sidx, ' binnary matrix file exist']);
    end
end

for k = 1:Nsub % save output
    subID = sublist{k};
    binconn_filepath = [outpath, '1_binconn/', subID];
    if ~exist([binconn_filepath,'_HO_binconn.mat'], 'file')
        FO_binconn = binconn_all{k}.FO_binconn;  % This pulls the structure for subject k
        Int_binconn = binconn_all{k}.Int_binconn;
        HO_binconn = binconn_all{k}.HO_binconn;

        % Save the structure in a file, typically using the subject ID in the filename
        save([binconn_filepath, '_FO_binconn.mat'], 'FO_binconn', '-v7.3');
        save([binconn_filepath, '_Int_binconn.mat'], 'Int_binconn', '-v7.3');
        save([binconn_filepath, '_HO_binconn.mat'], 'HO_binconn', '-v7.3');
        
        disp(['Subject ', subID, ' binary matrix successfully saved']);
    else
        disp(['Subject ', subID, ' binary matrix file exist']);
    end
end

%% w/o seed selection ver.
if ~exist([outpath,'1_binconn/sym_ver'],'dir') % creating folder
    mkdir([outpath,'1_binconn/sym_ver'])
end

for k = 1 : numel(sublist) 
    sublist_k = sublist{k};
    sidx = sublist_k;
    binconn_filepath = [outpath, '1_binconn/sym_ver/', sidx];
   
    if ~exist([binconn_filepath, '_binconn.mat'],'file')
        disp(['subject = ', num2str(sidx)])
        conn = load([outpath, '0_Connectivity_matrix/', sidx, '_WholeCTX_coef_FDR_corr.mat']).whl_thal_conn_matrix_correct;
        
        % binarize
        binconn = binarize_conn(conn);
        %binconn_all{k} = struct('sidx', sidx, 'binconn', binconn);
        save([binconn_filepath,'_binconn.mat'], 'binconn', '-v7.3');
        clear binconn
        disp(['Subject ', sidx, ' binnary matrix construction completed']);
    else
        disp(['Subject ', sidx, ' binnary matrix file exist']);
    end    
end

