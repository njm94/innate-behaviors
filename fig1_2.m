% Code to generate figure supplement 1-2. Num unilateral events for right 
% vs left sided stimulus. Also grooming over days

clear, clc
% This code will read the clustered behaviors and generate figure 1
addpath('/home/user/Documents/grooming/utils')
cluster_data = (fix_path('Y:\nick\behavior\grooming\20241114092737_behavior_clustering.mat'));

% map dendrogram labels to behaviors - Note this is subject to change based
% on the ordering of the output labels from the dendrogram in python
% script. Double check this.
label_map = ["Right", "Left", "Left Asymmetric", "Elliptical Left", ... 
    "Elliptical Right", "Right Asymmetric", "Elliptical"];

fileID = fopen('expt1_datalist.txt','r');

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = ''; 
% addpath('C:\Users\user\Documents\Nick\grooming\deprecated')

% Notes on 2p data
% ETR2 and ETR3 had corrupted videos for spontaneous trials, so we will
% just use evoked videos
%
% ECL3, IDR3, and RR3 have 3 days of spontaneous behavior - use last 2 days
% for the timeline plot. The spontaneous data may not have Boris labelling.
% Use the snippets data instead. Also these mice have longer baselines in
% the evoked trials. (5mins in these mice vs 30s in the 1p mice). 
% They don't do any grooming anyway so should be comparable
%
% ECL3 20240729 corrupted synchronization between brain and behavior. Can't
% use the brain data with grooming together but grooming alone should be ok

mp_list = {'Y:\nick\behavior\grooming\2p\ETR2_thy1\20231113143925'; 
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231113155903'; 
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231115174148'; 
    '/media/user/teamshare/nick/behavior/grooming/2p/ECL3_thy1/20240722/';
    '/media/user/teamshare/nick/behavior/grooming/2p/ECL3_thy1/20240724/';
    '/media/user/teamshare/nick/behavior/grooming/2p/ECL3_thy1/20240726/';
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240729'; 
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240731';
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240802';
    '/media/user/teamshare/nick/behavior/grooming/2p/IDR3_tTA6s/20240722/';
    '/media/user/teamshare/nick/behavior/grooming/2p/IDR3_tTA6s/20240724/';
    '/media/user/teamshare/nick/behavior/grooming/2p/IDR3_tTA6s/20240726/';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240729';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240731';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240802';
    '/media/user/teamshare/nick/behavior/grooming/2p/RR3_tTA8s/20240722/';
    '/media/user/teamshare/nick/behavior/grooming/2p/RR3_tTA8s/20240724/';
    '/media/user/teamshare/nick/behavior/grooming/2p/RR3_tTA8s/20240726/';
    'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240729';
    'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240802'};

data_list{1} = [data_list{1}; mp_list];
for i = 1:size(data_list{1},1)
    data_list{1}{i} = fix_path(data_list{1}{i});
end


% make behavior event raster
spon_index = [1,2,8,9,14,15,20,21,29:31,35:37,41:43];
evoked_index = ones(size(data_list{1}));
evoked_index(spon_index) = 0;
evoked_index = find(evoked_index);


%%
% get the drop info from boris file to determine if left vs right trial
trial_type = cell(1, length(data_list{1}));
include_boris = true; 
for i = 1:size(data_list{1},1)
    data_dir = data_list{1}{i};
    if any(i == spon_index)
        trial_type{i} = 'spon';
    else
        if ~isempty(getAllFiles(data_dir, '.tsv'))
            boris_file = [data_dir, filesep, getAllFiles(data_dir, '.tsv')];
            [events, snippets, b_table, video_end, cluster_labels] = get_labels_from_clustering_results(cluster_data, boris_file, include_boris);
    
            if isempty(events) % no grooming behaviors or drops
                error('No events')
            end
    
            if any(contains(events.Properties.VariableNames, 'DropRight'))
                num_drops_right = sum(events.DropRight);
            else
                num_drops_right = 0;
            end
    
            if any(contains(events.Properties.VariableNames, 'DropLeft'))
                num_drops_left = sum(events.DropLeft);
            else
                num_drops_left = 0;
            end
    
            if num_drops_right > num_drops_left
                trial_type{i} = 'right';
            elseif num_drops_left > num_drops_right
                trial_type{i} = 'left';
            elseif num_drops_right == num_drops_left
                disp(['Check this trial: ', num2str(i)])
                disp(['Manually determined to be right. Add drop labels in BORIS later'])
                trial_type{i} = 'right';
%                 continue
            end
    
        else % no boris file - use snippets
            disp('No Boris file')
            snippets_dir = [data_dir, filesep, 'snippets'];
            [snippets, ~] = parse_snippets(snippets_dir);

        end
    end
end

%% get index for right/left trials corresponding to cluster assignments
load(cluster_data);
bFiles = fix_path(cellstr(bFiles));

all_right = zeros(length(bFiles),1);
all_left = zeros(length(bFiles),1);

count_right = 1;
count_left = 1;
for i = 1:length(data_list{1})
    data_dir = fix_path(data_list{1}{i});

    matching_sessions = contains(bFiles, data_dir)';
    iTrial_type = trial_type{i};
    if strcmp(iTrial_type, 'right')
        all_right = all_right | matching_sessions;
        num_gR_dR(count_right) = numel(dendrogram_labels(matching_sessions & dendrogram_labels == find(strcmp(label_map, 'Right'))));
        num_gRA_dR(count_right) = numel(dendrogram_labels(matching_sessions & dendrogram_labels == find(strcmp(label_map, 'Right Asymmetric'))));
        num_gRE_dR(count_right) = numel(dendrogram_labels(matching_sessions & dendrogram_labels == find(strcmp(label_map, 'Elliptical Right'))));

        num_gL_dR(count_right) = numel(dendrogram_labels(matching_sessions & dendrogram_labels == find(strcmp(label_map, 'Left'))));
        num_gLA_dR(count_right) = numel(dendrogram_labels(matching_sessions & dendrogram_labels == find(strcmp(label_map, 'Left Asymmetric'))));
        num_gLE_dR(count_right) = numel(dendrogram_labels(matching_sessions & dendrogram_labels == find(strcmp(label_map, 'Elliptical Left'))));


        count_right = count_right + 1;
    elseif strcmp(iTrial_type, 'left')
        all_left = all_left | matching_sessions;

        num_gR_dL(count_left) = numel(dendrogram_labels(matching_sessions & dendrogram_labels == find(strcmp(label_map, 'Right'))));
        num_gRA_dL(count_left) = numel(dendrogram_labels(matching_sessions & dendrogram_labels == find(strcmp(label_map, 'Right Asymmetric'))));
        num_gRE_dL(count_left) = numel(dendrogram_labels(matching_sessions & dendrogram_labels == find(strcmp(label_map, 'Elliptical Right'))));

        num_gL_dL(count_left) = numel(dendrogram_labels(matching_sessions & dendrogram_labels == find(strcmp(label_map, 'Left'))));
        num_gLA_dL(count_left) = numel(dendrogram_labels(matching_sessions & dendrogram_labels == find(strcmp(label_map, 'Left Asymmetric'))));
        num_gLE_dL(count_left) = numel(dendrogram_labels(matching_sessions & dendrogram_labels == find(strcmp(label_map, 'Elliptical Left'))));

        count_left = count_left + 1;
    end
end

%%
close all


figure
boxplot([num_gR_dR' num_gL_dR'], 'colors', 'k', 'symbol', '')
hold on
hline(0, 'k--')
swarmchart([ones(size(num_gR_dR))-0.25 1.25+ones(size(num_gR_dR))], [num_gR_dR num_gL_dR], [], 'k', 'XJitterWidth', 0.25)
if isnormal(num_gR_dR) && isnormal(num_gL_dR)
    [~,p] = ttest(num_gR_dR, num_gL_dR);
else
    disp('not normal')
    p = ranksum(num_gR_dR, num_gL_dR);
end

title(num2str(p))
ylabel('Right stimulus')
xticklabels({'Right', 'Left'})
xlabel([num2str(median(num_gR_dR)), '+/-', num2str(iqr(num_gR_dR)), '-------',num2str(median(num_gL_dR)), '+/-', num2str(iqr(num_gL_dR))])
axis([.5 2.5 -5 150])
%%
figure
boxplot([num_gRA_dR' num_gLA_dR'], 'colors', 'k', 'symbol', '')
hold on
hline(0, 'k--')
swarmchart([ones(size(num_gRA_dR))-0.25 1.25+ones(size(num_gRA_dR))], [num_gRA_dR num_gLA_dR], [], 'k', 'XJitterWidth', 0.25)
if isnormal(num_gRA_dR) && isnormal(num_gLA_dR)
    [~,p] = ttest(num_gRA_dR, num_gLA_dR);
else
    disp('Not normal')
    p = ranksum(num_gRA_dR, num_gLA_dR);
end
title(num2str(p))
xticklabels({'Right Asymmetric', 'Left Asymmetric'})
xlabel([num2str(median(num_gRA_dR)), '+/-', num2str(iqr(num_gRA_dR)), '-------',num2str(median(num_gLA_dR)), '+/-', num2str(std(num_gLA_dR))])
ylabel('Right stimulus')
axis([.5 2.5 -5 150])
%%
figure
boxplot([num_gRE_dR' num_gLE_dR'], 'colors', 'k', 'symbol', '')
hold on
hline(0, 'k--')
swarmchart([ones(size(num_gRE_dR))-0.25 1.25+ones(size(num_gRE_dR))], [num_gRE_dR num_gLE_dR], [], 'k', 'XJitterWidth', 0.25)
if isnormal(num_gRE_dR) && isnormal(num_gLE_dL)
    [~,p] = ttest(num_gRE_dR, num_gLE_dR);
else
    disp('Not normal')
    p = ranksum(num_gRE_dR, num_gLE_dR);
end
title(num2str(p))
xticklabels({'Elliptical Right', 'Elliptical Left'})
xlabel([num2str(median(num_gRE_dR)), '+/-', num2str(std(num_gRE_dR)), '-------',num2str(median(num_gLE_dR)), '+/-', num2str(std(num_gLE_dR))])
ylabel('Right stimulus')
axis([.5 2.5 -5 150])
%%
figure
boxplot([num_gR_dL' num_gL_dL'], 'colors', 'k', 'symbol', '')
hold on
swarmchart([ones(size(num_gR_dL))-0.25 1.25+ones(size(num_gR_dL))], [num_gR_dL num_gL_dL], [], 'k', 'XJitterWidth', 0.25)
hline(0, 'k--')
if isnormal(num_gR_dL) && isnormal(num_gL_dL)
    [~,p] = ttest(num_gR_dL, num_gL_dL);
else
    disp('not normal')
    p = ranksum(num_gR_dL, num_gL_dL);
end
title(num2str(p))
xticklabels({'Right', 'Left'})
xlabel([num2str(median(num_gR_dL)), '+/-', num2str(std(num_gR_dL)), '-------',num2str(median(num_gL_dL)), '+/-', num2str(std(num_gL_dL))])
ylabel('Left stimulus')
axis([.5 2.5 -5 150])
%%
figure
boxplot([num_gRA_dL' num_gLA_dL'], 'colors', 'k', 'symbol', '')
hold on
swarmchart([ones(size(num_gRA_dL))-0.25 1.25+ones(size(num_gRA_dL))], [num_gRA_dL num_gLA_dL], [], 'k', 'XJitterWidth', 0.25)
if isnormal(num_gRA_dL) && isnormal(num_gLA_dL)
    [~,p] = ttest(num_gRA_dL, num_gLA_dL);
else
    disp('Not normal')
    p = ranksum(num_gRA_dL, num_gLA_dL);
end
hline(0, 'k--')
title(num2str(p))
xticklabels({'Right Asymmetric', 'Left Asymmetric'})
xlabel([num2str(median(num_gRA_dL)), '+/-', num2str(std(num_gRA_dL)) '-------',num2str(median(num_gLA_dL)), '+/-', num2str(std(num_gLA_dL))])
ylabel('Left stimulus')
axis([.5 2.5 -5 150])
%%
figure
boxplot([num_gRE_dL' num_gLE_dL'], 'colors', 'k', 'symbol', '')
hold on
hline(0, 'k--')
swarmchart([ones(size(num_gRE_dL))-0.25 1.25+ones(size(num_gRE_dL))], [num_gRE_dL num_gLE_dL], [], 'k', 'XJitterWidth', 0.25)
if isnormal(num_gRE_dL) && isnormal(num_gLE_dL)
    [~,p] = ttest(num_gRE_dL, num_gLE_dL);
else
    disp('Not normal')
    p = ranksum(num_gRE_dL, num_gLE_dL);
end
title(num2str(p))
xticklabels({'Elliptical Right', 'Elliptical Left'})
xlabel([num2str(median(num_gRE_dL)),'+/-', num2str(std(num_gRE_dL)) '-------',num2str(median(num_gLE_dL)), '+/-', num2str(std(num_gLE_dL))])
ylabel('Left stimulus')
axis([.5 2.5 -5 150])

%%
data_plots = {'DropRightUnilateral', 'DropRightAsymmetric', 'DropRightEllipAsymm', 'DropLeftUnilateral', 'DropLeftAsymmetric', 'DropLeftEllipAsymm'};
for i = 1:6
    ax = figure(i);
    saveas(ax, fix_path(['Y:\nick\behavior\grooming\figures\',data_plots{i}, '.svg']))
end