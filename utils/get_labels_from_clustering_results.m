function [events, b_idx, t, video_end, cluster_labels] = get_labels_from_clustering_results(cluster_data, boris_tsv, include_boris)
% Get labels from UMAP clustering analysis in the same format as read_boris
%
% Inputs: 
%   cluster_data    (full path to UMAP cluster results)
%   boris_file      (current mouse BORIS file. Still use this for lick and
%                   video end information)
%   include_boris   (include lick and drop information from BORIS file into
%                   event matrix. Default FALSE)
%
% Outputs:
%       events          (Binary table [trial_length x num_behaviors])
%       b_idx           (Cell array of event indices [trial_length x 1])
%       t               (Boris raw table)
%       video_end       (Index of video end from BORIS label)
%       cluster_labels  (Hierarchical cluster label)
%
% Usage: 
%       [events, b_idx, t] = read_boris(boris_tsv);
%       [events, b_idx, t] = read_boris(boris_tsv, len);

if nargin<3 || isempty(include_boris), include_boris = false; end

% load umap clustering results
% umap_file = fix_path('Y:\nick\behavior\grooming\220241108104909_behavior_clustering.mat');
umap_results = load(cluster_data);
labels = umap_results.dendrogram_labels;
b_idx = umap_results.bIdx;

% map dendrogram labels to behaviors - Note this is subject to change based
% on the ordering of the output labels from the dendrogram in python
% script. Double check this.
% label_map = ["Right", "Left", "Left Asymmetric", ...
%     "Elliptical Asymmetric", "Right Asymmetric", "Elliptical"];
label_map = ["Right", "Left", "Left Asymmetric", "Elliptical Left", ...
    "Elliptical Right", "Right Asymmetric", "Elliptical"];

% get video end information from boris
[boris_events, ~, t, video_end] = read_boris(boris_tsv);


% Match behavior files using the boris filename
for i = 1:size(umap_results.bFiles,1)
    if isunix
        newBfiles(i,:) = strrep(umap_results.bFiles(i,:), 'teamshare/nick', 'teamshare/TM_Lab/nick');
        % newBfiles = umap_results.bFiles;
    else
        if contains(umap_results.bFiles(i,:), '/media/user/teamshare/TM_Lab/')
            newBfiles(i,:) = strrep(umap_results.bFiles(i,:), '/media/user/teamshare/TM_Lab/', 'Y:\');
        elseif contains(umap_results.bFiles(i,:), '/media/user/teamshare/')
            newBfiles(i,:) = strrep(umap_results.bFiles(i,:), '/media/user/teamshare/', 'Y:\');
        else
            disp(['Check filepaths for ', umap_results.bFiles(i,:)])
        end
        newBfiles(i,:) = strrep(newBfiles(i,:), '/', '\');
    end
end


bFiles = cellstr(newBfiles);
clear newBfiles

matchingIndex = strcmp(bFiles, boris_tsv);

% get the clustering labels and behavior index for the correct matched
% video
labels = labels(matchingIndex);
b_idx = b_idx(matchingIndex,:);


% create the events table
events = zeros(video_end, length(label_map));
for i = 1:size(b_idx,1)
    events(b_idx(i,1):b_idx(i,2),strcmp(label_map, label_map(labels(i)))) = 1;
end
empty_idx = sum(events)==0;
events = array2table(events, 'VariableNames', label_map);

cluster_labels = label_map(labels);


% remove empty variables
events = removevars(events, empty_idx);

% add lick and drop events to event matrix if specified
if include_boris 
    % consolidate lick events into a single variable in the table
    lick_idx = contains(boris_events.Properties.VariableNames, 'Lick');
    Lick = boris_events(:,lick_idx);
    Lick = any(table2array(Lick),2); 
    if any(Lick)
        events = addvars(events, Lick);
    end

    drop_idx = find(contains(boris_events.Properties.VariableNames, 'Drop'));
    for i = 1:length(drop_idx)
        if contains(boris_events.Properties.VariableNames(drop_idx(i)), 'Left')
            DropLeft = table2array(boris_events(:,drop_idx(i)));
            events = addvars(events, DropLeft);
        elseif contains(boris_events.Properties.VariableNames(drop_idx(i)), 'Right')
            DropRight = table2array(boris_events(:,drop_idx(i)));
            events = addvars(events, DropRight);
        else
            DropCenter = table2array(boris_events(:,drop_idx(i)));
            events = addvars(events, DropCenter);
        end  
    end
  
end


end