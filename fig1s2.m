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
include_boris = true; 
for i = 5%1:size(data_list{1},1)
    data_dir = data_list{1}{i};
    if ~isempty(getAllFiles(data_dir, '.tsv'))
        boris_file = [data_dir, filesep, getAllFiles(data_dir, '.tsv')];
        [events, snippets, b_table, video_end, cluster_labels] = get_labels_from_clustering_results(cluster_data, boris_file, include_boris);
        gadadg
        if isempty(events) % no grooming behaviors
            continue
        end

        if contains(events.Properties.VariableNames, 'DropRight')
            num_drops_right = sum(events.DropRight);
        else
            num_drops_right = 0;
        end

        if contains(events.Properties.VariableNames, 'DropLeft')
            num_drops_left = sum(events.DropLeft);
        else
            num_drops_left = 0;
        end

        if num_drops_right > num_drops_left
            trial_type = 'right';
        elseif num_drops_left > num_drops_right
            trial_type = 'left';
        elseif num_drops_right == num_drops_left
            if num_drops_right + num_drops_left == 0
                trial_type = 'spon';
            else
                error('Check this trial')
            end
        end

    else % no boris file - use snippets
        disp('No Boris file')
        [snippets, labels] = parse_snippets(snippets_dir);
    end
end
