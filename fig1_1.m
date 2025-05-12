clear
load(fix_path('Y:\nick\behavior\grooming\20241114092737_behavior_clustering.mat'))

figure, hold on

label_map = {'Right', 'Left', 'Left Asymmetric', 'Elliptical Left', 'Elliptical Right', 'Right Asymmetric', 'Elliptical'};
labs = unique(cellstr(manual_labels));
manual_label_map = {'Elliptical', 'Elliptical Asymmetric', 'Large Bilateral', 'Left', 'Left Asymmetric', 'Right', 'Right Asymmetric'};
for i = 1:length(labs)
%     scatter(embedding(dendrogram_labels==labs(i),1), embedding(dendrogram_labels==labs(i),2),'.')
    scatter(embedding(strcmp(cellstr(manual_labels), labs{i}), 1), embedding(strcmp(cellstr(manual_labels), labs{i}), 2), '.')
end
legend(manual_label_map, 'Box','off')
%%
ax = gca;
saveas(ax, fix_path(['Y:\nick\behavior\grooming\figures\','manual_UMAP', '.svg']))


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
clc
behavior_files = unique(bFiles);
labs = unique(dendrogram_labels);


counts = cell(1, length(labs));
% figure

xEdges = linspace(-200, 200, 50);
yEdges = linspace(-220, 180, 50);

for i = 1:size(behavior_files,1)
    % load DLC tracks
    data_root = fix_path(fileparts(behavior_files{i}));
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
        if isempty(counts{j}), counts{j} = zeros(length(xEdges)-1); end
        figure(j), hold on
        for k = 1:length(clus_j)
            plot(-flr_x(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), flr_y(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), 'Color', [1 0 0 0.1])
            plot(-fll_x(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), fll_y(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), 'Color', [0 0 1 0.1])
            counts{j} = counts{j} + histcounts2(flr_y(bIdx(clus_j(k),1):bIdx(clus_j(k),2)),-flr_x(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), yEdges, xEdges);
            counts{j} = counts{j} - histcounts2(fll_y(bIdx(clus_j(k),1):bIdx(clus_j(k),2)),-fll_x(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), yEdges, xEdges);
        end
        axis([-200 200 -220 180])
        xticks([])
        yticks([])
        box on
    end
end
%%
for i = 1:length(labs)
    exportgraphics(figure(i), fix_path(append('Y:\nick\behavior\grooming\figures\', label_map(i),'.png')), 'Resolution', 300, 'ContentType', 'image')
end
% fig = gcf;
% print(fig,fix_path(['Y:\nick\behavior\grooming\figures\','alltracks', '.svg']),'-dsvg');

%%


for i = 1:length(counts)
    figure
    tmp = counts{i};
    tmp(tmp>0) = tmp(tmp>0)./max(tmp(:));
    tmp(tmp<0) = -tmp(tmp<0)./min(tmp(:));
    imagesc(imgaussfilt(tmp, 1))
    set(gca,'YDir','normal')
    caxis([-1 1])
    colormap(bluewhitered())
    xlim([-6.0634 56.0634])
    ylim([0.5000 49.5000])
    xticks([])
    yticks([])
    
    
    ax = gca;
%     exportgraphics(ax, fix_path(append('Y:\nick\behavior\grooming\figures\', label_map(i),'.png')), 'Resolution', 300, 'ContentType', 'image')
end


