figure, hold on
labs = unique(dendrogram_labels);
for i = 1:length(labs)
    scatter(embedding(dendrogram_labels==labs(i),1), embedding(dendrogram_labels==labs(i),2),'.')
end


%%
for i = 1:size(bFiles,1)
    if isunix
        newBfiles = bFiles;
    else
        newBfiles(i,:) = strrep(bFiles(i,:), '/media/user/teamshare', 'Y:');
        newBfiles(i,:) = strrep(newBfiles(i,:), '/', '\');
    end
end

bFiles = cellstr(newBfiles);
%%

behavior_files = unique(bFiles);
figure
for i = 1:size(behavior_files,1)
    % load DLC tracks
    data_root = fileparts(behavior_files(i,:));
    dlc_pos = readmatrix([data_root, filesep, getAllFiles(data_root, '1030000.csv')]);
    nose_x = median(dlc_pos(:,1));
    nose_y = median(dlc_pos(:,2));
    
    % get all position data relative to nose to account for variability
    % in camera positioning across trials
    flr_x = nose_x - dlc_pos(:,4);
    flr_y = nose_y - dlc_pos(:,5);    
    fll_x = nose_x - dlc_pos(:,7);
    fll_y = nose_y - dlc_pos(:,8);
    
    matching_sessions = contains(bFiles, data_root);
    snippets_dir = [data_root filesep 'snippets'];
    for j = 1:length(labs)
        clus_j = find(dendrogram_labels==labs(j) & matching_sessions');
        if isempty(clus_j), continue; end
        subplot(1,length(labs), j), hold on
        for k = 1:length(clus_j)
            plot(-flr_x(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), flr_y(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), 'Color', [1 0 1 0.1])
            plot(-fll_x(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), fll_y(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), 'Color', [0 1 1 0.1])
        end
    end
end