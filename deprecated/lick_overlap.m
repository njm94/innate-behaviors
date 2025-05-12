%% Lick overlap analysis with other behaviors



clc, clear, close all
fileID = fopen('expt1_datalist.txt','r');

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = ''; 

mp_list = {'Y:\nick\behavior\grooming\2p\ETR2_thy1\20231113143925'; ...
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231113155903'; ...
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231115174148'; ...
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240729'; ...
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240731';
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240802';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240729';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240731';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240802';
    'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240729';
    'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240802'};

cluster_data = fix_path('Y:\nick\behavior\grooming\20241114092737_behavior_clustering.mat');

data_list = [data_list{1}; mp_list];

%%

include_boris = true;
N = length(data_list);

all_behaviors = {'Right', 'Left', 'Elliptical', 'Right Asymmetric', ...
    'Left Asymmetric', 'Elliptical Right', 'Elliptical Left'};
lick_overlap_counter = zeros(1, length(all_behaviors)+1);
groom_overlap_counter = zeros(1, length(all_behaviors));
groom_totals = zeros(1, length(all_behaviors));

buffer_size = 0; % frame buffer to add before/after each lick or groom event for determining the overlap. this is 1/6 of a second
    for j = 1:N
        data_dir = fix_path(data_list{j});
        timestamp_file = [data_dir, filesep, getAllFiles(data_dir, '_trim.txt')];
        [mouse_root_dir, exp_date, ~] = fileparts(data_dir);
        [~, mouse_id, ~] = fileparts(mouse_root_dir);
        [expt_root_dir, ~, ~] = fileparts(mouse_root_dir);
    
    
        boris_file = [data_dir, filesep, getAllFiles(data_dir, '.tsv')];
        [events, snippets, b_table, video_end, cluster_labels] = get_labels_from_clustering_results(cluster_data, boris_file, include_boris);
        if isempty(events) 
            disp('Events matrix is empty. skipping...')
            continue
        end
        labels = events.Properties.VariableNames;
    
        if ~any(contains(labels, 'Lick'))
            disp('No licks detected. Skipping...')
            continue
        end
        
        lick = events.Lick;
        lick_idx = arr2idx(lick);
    
        idx = contains(events.Properties.VariableNames, 'Lick') | ...
        contains(events.Properties.VariableNames, 'Drop') | ...
        contains(events.Properties.VariableNames, 'Video') | ...
        contains(events.Properties.VariableNames, 'Flail') ;
        groom_events = removevars(events, idx);

        %% Get the total counts for grooming events
        % This needs to be a separate loop to avoid double counting groom
        % events

        for ii = 1:size(snippets,1)
            current_behavior = cluster_labels(ii);
            behavior_index = find(strcmp(all_behaviors, current_behavior));
            if any(lick(snippets(ii,1)-buffer_size:min([snippets(ii,2)+buffer_size, length(lick)])))
                groom_overlap_counter(1, behavior_index) = groom_overlap_counter(1, behavior_index) + 1;
            end           
            groom_totals(1, behavior_index) = groom_totals(1, behavior_index) + 1;
        end

        %% Find the proportion of licks that overlap with each grooming behavior
        % Iterate through lick events to count grooming behaviors that overlap
        % with licks Store overlap counts in an array of size [N+1, 1] where N
        % is the number of grooming behaviors (+1 for the count of 
        % non-overlapping licks)    
        lick_already_counted = false;
        for i = 1:size(lick_idx,1)
            for ii = 1:size(groom_events,2)
                tmp = table2array(groom_events(:,ii));
                current_behavior = events.Properties.VariableNames{ii};
                behavior_index = find(strcmp(all_behaviors, current_behavior));
      
                % if there is an overlap detected, count it and then break out
                % of the the loop. Since all grooming behaviors are 
                % non-overlapping, there is no point in iterating
                % through the other behaviors.
                % min max statements are to account for licks that occur
                % near the start/end of the session
                if any(tmp(max([lick_idx(i,1)-buffer_size, 1]):min([lick_idx(i,2)+buffer_size, length(lick)])))
                    lick_overlap_counter(1, behavior_index) = lick_overlap_counter(1, behavior_index) + 1;
                    lick_already_counted = true;
                    break
                end            
            end
    
            if lick_already_counted
                % reset
                lick_already_counted = false;
            else
                lick_overlap_counter(1, end) = lick_overlap_counter(1, end) + 1;
            end
        end
    end


%%
figure, bar(100*lick_overlap_counter./sum(lick_overlap_counter,2)), 
xticklabels([all_behaviors, {'None'}])
ylabel('Percentage')
title('Licks overlapping with grooming behaviors')

figure, bar(100*groom_overlap_counter./groom_totals), 
xticklabels(all_behaviors)
ylabel('Percentage')
title('Grooming behaviors overlapping with licks')
%%




