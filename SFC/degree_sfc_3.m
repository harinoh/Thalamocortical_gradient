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
Nstep = 200;

%% 3) Save DC values per ROI/Network for all subjects
seed_name = 'HO_binconn';
seed = strsplit(seed_name,'_');
seed = seed{1};

%net = MRIread([basepath,'Gradient_result/3_2_FunctionalNetwork/Thalamic_FN_Gordon_Bi_thal_222_SCZOHC_FWHM4.nii.gz']).vol;
%num_network = 7;

roi_dc = zeros(Nsub, Nroi, Nstep);
net_dc = zeros(Nsub, Nstep);
for sidx = 1 : Nsub
    subID = sublist{sidx};
    disp(['Loading subject ', subID, ' for analysis']);
    sfc = load([outpath, '2_sfc/',seed_name,'/', subID,'_sfc.mat']).sfc;
    for step = 1 : Nstep
        dc = sum(sfc(:,:,step), 1);
        dc(isinf(dc)|isnan(dc)) = 0;
        roi_dc(sidx, :, step) = dc;
        net_dc(sidx, step) = mean(dc);
        %for nidx = 1 : num_network
        %    if any(net == nidx)
        %        net_dc(sidx, nidx, step) = mean(dc(net == nidx));
        %    end
        %end
    end
    disp(['Subject ', subID, ' degree of centrality calculated']);
    clear sfc
end
save([outpath, '2_sfc/wholesub_',seed,'_ROI_dc.mat'], 'roi_dc', '-v7.3');
save([outpath, '2_sfc/wholesub_',seed,'_NET_dc.mat'], 'net_dc', '-v7.3');
end
