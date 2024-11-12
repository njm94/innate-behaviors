function [events, b_idx, t, video_end, cluster_labels] = get_labels_from_clustering_results(cluster_data, boris_tsv)
% Get labels from UMAP clustering analysis in the same format as read_boris
%
% Inputs: 
%   cluster_data    (full path to UMAP cluster results)
%   boris_file      (current mouse BORIS file. Still use this for lick and
%                   video end information)
%
% Outputs:
%       events          (Binary table [trial_length x num_behaviors])
%       b_idx           (Cell array of event indices [trial_length x 1])
%       t               (Boris raw table)
%
% Usage: 
%       [events, b_idx, t] = read_boris(boris_tsv);
%       [events, b_idx, t] = read_boris(boris_tsv, len);


% load umap clustering results
% umap_file = fix_path('Y:\nick\behavior\grooming\220241108104909_behavior_clustering.mat');
umap_results = load(cluster_data);

% label map
label_map = ["Right", "Left", "Elliptical Asymmetric", ...
    "Left Asymmetric", "Right Asymmetric", "Elliptical"];

% get video end information from boris
[~, ~, t, video_end] = read_boris(boris_tsv);


% Match behavior files using the boris filename
for i = 1:size(umap_results.bFiles,1)
    if isunix
        newBfiles = umap_results.bFiles;
    else
        newBfiles(i,:) = strrep(umap_results.bFiles(i,:), '/media/user/teamshare', 'Y:');
        newBfiles(i,:) = strrep(newBfiles(i,:), '/', '\');
    end
end
bFiles = cellstr(newBfiles);
matchingIndex = strcmp(bFiles, boris_tsv);

% get the clustering labels and behavior index for the correct matched
% video
labels = umap_results.dendrogram_labels(matchingIndex);
b_idx = umap_results.bIdx(matchingIndex,:);

% create the events table
events = zeros(video_end, length(label_map));
for i = 1:size(b_idx,1)
    events(b_idx(i,1):b_idx(i,2),strcmp(label_map, label_map(labels(i)))) = 1;
end
empty_idx = sum(events)==0;
events = array2table(events, 'VariableNames', label_map);

cluster_labels = label_map(umap_results.dendrogram_labels(matchingIndex));


% remove empty variables
events = removevars(events, empty_idx);

end