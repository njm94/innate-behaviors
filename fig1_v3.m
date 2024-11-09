% This code will read the clustered behaviors and generate figure 1

load('Y:\nick\behavior\grooming\20241108104909_behavior_clustering.mat')

figure, hold on
% show UMAP embedding

labs = unique(dendrogram_labels);
for i = 1:length(labs)
    scatter(embedding(dendrogram_labels==labs(i),1), embedding(dendrogram_labels==labs(i),2),'.')
    axis off
end

%%
% map dendrogram labels to behaviors
label_map = ["Right", "Left", "Elliptical Asymmetric", ...
    "Left Asymmetric", "Right Asymmetric", "Elliptical"];

clc, close all
% get the example images
example_img_path = fix_path('Y:\nick\behavior\grooming\grooming_example_images');
example_images = getAllFiles(example_img_path);
data_root = fix_path('Y:\nick\behavior\grooming\1p\ECR2_thy1');
for i = 1:length(label_map)
    label_str = append(lower(strrep(label_map(i), ' ', '_')), '.');
    img_path = example_images{contains(example_images, label_str)};
    img = loadtiff([example_img_path, filesep, img_path]);

    % get nose for that particular file as reference
    tmp = strfind(img_path, '_');
    exp_date = img_path(tmp(1)+1:tmp(2)-1);
    dlcpos_file = getAllFiles([data_root, filesep, exp_date], '1030000.csv');
    dlc_pos = readmatrix([data_root, filesep, exp_date, filesep, dlcpos_file]);

    nose_ref_x(i) = median(dlc_pos(:,1));
    nose_ref_y(i) = median(dlc_pos(:,2));
    
    % Brighten image so it's easier to see.
    % All images taken from one session except Elliptical Asymmetric, which
    % didn't occur in that session. The EA behavior image is slightly
    % brighter than the others, so adjustment is a bit different
    img = log(double(img(nose_ref_y(i)-100:nose_ref_y(i)+200, nose_ref_x(i)-150:nose_ref_x(i) + 150)));
    img = 255*img./max(img(:));
    if strcmpi(label_map(i), 'Elliptical Asymmetric')
        img(img<130) = 130;
    else
        img(img<100) = 100;
        img(img>225) = 225;
    end
    img = 255*(img-min(img(:)))./(max(img(:))-min(img(:)));
    figure, imagesc(img), hold on, colormap gray
    
    axis equal off

end

% Deal with file paths
for i = 1:size(bFiles,1)
    if isunix
        newBfiles = bFiles;
    else
        newBfiles(i,:) = fix_path(bFiles(i,:));
    end
end
bFiles = cellstr(newBfiles);
%%

behavior_files = unique(bFiles);

counter = zeros(length(labs),1);

for i = 1:size(behavior_files,1)
    % load DLC tracks
    data_root = fileparts(behavior_files(i,:));
    dlc_pos = readmatrix([data_root, filesep, getAllFiles(data_root, '1030000.csv')]);
    nose_x = median(dlc_pos(:,1));
    nose_y = median(dlc_pos(:,2));
    
    % First, get all position data relative to nose to account for 
    % variability in camera positioning across trials
    % Then, use reference nose position for plotting in the final images
    flr_x = nose_x - dlc_pos(:,4);
    flr_y = dlc_pos(:,5) - nose_y;    
    fll_x = nose_x - dlc_pos(:,7);
    fll_y = dlc_pos(:,8) - nose_y;
    
    matching_sessions = contains(bFiles, data_root);
    snippets_dir = [data_root filesep 'snippets'];
    for j = 1:length(labs)
        clus_j = find(dendrogram_labels==labs(j) & matching_sessions');
        if isempty(clus_j), continue; end
        
        figure(j),
        for k = 1:length(clus_j)
            counter(j) = counter(j)+1;
            if counter(j)>200, continue; end
            plot(nose_ref_x(j)-151-flr_x(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), nose_ref_y(j)+flr_y(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), 'Color', [1 0 1 0.15])
            plot(nose_ref_x(j)-151-fll_x(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), nose_ref_y(j)+fll_y(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), 'Color', [0 1 1 0.15])
        end
    end
end


%%

for i = 1:length(label_map)
    figure(i);
    axis([0, 300, 0, 300])
    ax = gca;
%     exportgraphics(ax, append('Y:\nick\behavior\grooming\figures', label_map(i),'.eps'), 'ContentType', 'vector')
    exportgraphics(ax, append('Y:\nick\behavior\grooming\figures\', label_map(i),'.png'), 'Resolution', 300)

end



