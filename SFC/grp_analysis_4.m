%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Perform group analysis using degree of SFC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grp_analysis_4_4()
clear all;
clc;

basepath = '/Volume/CCNC/harin_oh/1_thalamocortical/';
outpath = '/Volume/CCNC/harin_oh/1_thalamocortical/SFC_result/';
group_name = 'SCZOHC';
sublist = importdata([basepath, 'SCZOHC/Subject_list_',group_name,'.txt']); %change group

Nroi = 56447;
Nstep = 200;

%% 0) Divide participants into two group based on obesity phenotypes
demo = readtable([basepath, 'SCZOHC_rsfMRI_Demo.xlsx'],'Sheet','baseline_n129_4SFC','Range','A1:I130');
group2 = demo{:,6} +1; 
group = [demo{:,4},num2cell(group2)]; % Norm first
clear group2 demo

%% 1) Compute group average SFC
seed_name = 'Int';
% seed_idx = load([basepath, 'SFC_result/ctxthal_', seed_name, '_seed_region.mat']).comb_seed_idx_int;
group_ids = cell2mat(group(:,2));

grpmean_SFC = cell(2,1);
grpmean_SFC{1} = zeros(657, Nroi, Nstep); % grp1
grpmean_SFC{2} = zeros(657, Nroi, Nstep); % grp2
clear seed_idx

for k = 1 : numel(sublist)
    sidx = sublist{k};
    disp(['subject = ', sidx])  
    % Iterate over each row in X
    i = find(strcmp(group,sidx));
    if ~isempty(i)
%         group_idx = group{i,2}; % If there's a match, extract the value from the second column
        sfc = load([outpath, '2_sfc/',seed_name,'_binconn/',sidx,'_sfc.mat']).sfc;
        % divide the sfc into half
        sfc_1 = sfc(701:1357,:,:);
        clear sfc
        % Update the correct group
        grpmean_SFC{1} = grpmean_SFC{1} + sfc_1;
        clear sfc_1
    end
end
% Reload the final added SFC matrix
disp(['Calculating mean'])
grpmean_SFC{1} = grpmean_SFC{1} / sum(group_ids == 1);
grpmean_SFC{2} = grpmean_SFC{2} / sum(group_ids == 2);

% Calculating to degree of centrality
grpmean_DC = zeros(2, Nroi, Nstep);
grpmean_DC(1,:,:) = squeeze(sum(grpmean_SFC{1}(:,:,1:Nstep), 1));
grpmean_DC(2,:,:) = squeeze(sum(grpmean_SFC{2}(:,:,1:Nstep), 1));

save([outpath, '2_sfc/groupmean_SFC_',seed_name,'_part2.mat'], 'grpmean_SFC', '-v7.3');
save([outpath, '2_sfc/groupmean_SFC_DC_',seed_name,'_part2.mat'],'grpmean_DC', '-v7.3');
end
