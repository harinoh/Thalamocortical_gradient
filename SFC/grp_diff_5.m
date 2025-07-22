%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Perform group difference using degree of SFC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grp_diff_4_5()
outpath = '/Volume/CCNC/harin_oh/1_thalamocortical/SFC_result/';

Nroi = 56447;
Nstep = 200;
seed = 'FO'; % 'FO' 'Int' 'HO'

if ~exist([outpath,'3_stats/'],'dir') % creating folder
    mkdir([outpath,'3_stats/'])
end
%% 4) Group difference test: roi-level

roi_dc = load([outpath, '2_sfc/wholesub_',seed,'_ROI_dc.mat']).roi_dc;
H = zeros(Nroi, Nstep);
P = zeros(Nroi, Nstep);
T = zeros(Nroi, Nstep);

for step = 1 : Nstep
    grp1_dc = roi_dc(1:69,:,step);
    grp2_dc = roi_dc(70:129,:,step);
    for roi = 1 : Nroi
        [~,p,~,stats] = ttest2(grp1_dc(:,roi), grp2_dc(:,roi));
        P(roi, step) = p;
        T(roi, step) = stats.tstat;
    end
    % correct across networks
    [selected, ~, ~, corrected_P] = fdr_bh(P(:, step), 0.05);
    H(:, step) = selected;
    C_P(:, step) = corrected_P;
end
    
[FDR_index, FDRthd] = FDR(P(:,step), 0.05);
FDR_idx(:, step) = FDR_index;


save([outpath, '3_stats/groupdiff_ROI_',seed,'ttest.mat'], 'H', 'P', 'T', 'C_P', 'FDR_idx');

%% 5) Group difference test: network-level

net_dc = load([outpath, '2_sfc/wholesub_',seed,'_NET_dc.mat']).net_dc;
H = zeros(Nroi, Nstep);
P = zeros(Nroi, Nstep);
T = zeros(Nroi, Nstep);

for step = 1 : Nstep
    grp1_net_dc = net_dc(1:69,step);
    grp2_net_dc = net_dc(70:129,step);
    [~,p,~,stats] = ttest2(grp1_net_dc, grp2_net_dc);
    NET_P(step, 1) = p;
    NET_T(step, 1) = stats.tstat;
    
    % correct across networks
%     [selected, ~, ~, corrected_P] = fdr_bh(NET_P(:, step), 0.05);
%     NET_H(:, step) = selected;
%     NET_C_P(:, step) = corrected_P;
end
    
save([outpath, '3_stats/groupdiff_NET_',seed,'ttest.mat'], 'NET_H', 'NET_P', 'NET_T', 'NET_C_P');

end
