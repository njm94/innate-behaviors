% get avg behavior frames into form for UMAP clustering

clear, clc, close all
addpath(fix_path('C:\Users\user\Documents\Nick\grooming\utils'))
data_root = fix_path('Y:\nick\behavior\grooming\1p');
mice = {'ECR2_thy1', 'GER2_ai94', 'HYL3_tTA', 'IBL2_tTA'};

bvars_for_UMAP = ["Right", "Left", "Elliptical", "Right Asymmetric", ...
    "Left Asymmetric", "Elliptical Right", "Elliptical Left"];%, "Lick"];
load('allenDorsalMap.mat');
data_for_UMAP = [];
labels_for_UMAP = [];
for j = 1:length(mice)
    % load mask and transformation matrix to get data in Allen atlas space
    % also transform the mask into Allen atlas space
    load(fix_path([data_root, filesep, mice{j}, filesep, 'mask.mat']))
    load(fix_path([data_root, filesep, mice{j}, filesep, 'atlas_tform.mat']));
    % mask = imwarp(mask, tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));

    dff_path = [data_root, filesep, mice{j}, filesep, 'outputs'];
    dff_mat = getAllFiles(dff_path, 'avg_behavior_frames.mat');
    % If there are multiple versions, use the most recent one
    if size(dff_mat,1) > 1
        dff_mat = sort(dff_mat);
        dff_mat = dff_mat{end};
    end
    load([dff_path, filesep, dff_mat], 'avg_behavior_frames', 'bvars');

    for i = 1:length(bvars_for_UMAP)
        cell_idx = find(strcmp(bvars, bvars_for_UMAP(i)));

        % transform the behavior frames so they are in the reference space
        % of the Allen atlas
        tmp = mask.*(mean(pagetranspose(avg_behavior_frames{cell_idx}),3));
        % tmp = mask.*imwarp(tmp, tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
        data_for_UMAP = cat(1, data_for_UMAP, (reshape(tmp, size(tmp,1)*size(tmp,2), [])'));
        labels_for_UMAP = cat(1, labels_for_UMAP, cellstr(repmat(bvars_for_UMAP(i), size(tmp,3), 1)));


    end

end


%%
[coeff, score, ~, ~, explained] = pca(data_for_UMAP);
disp('Done with PCA')
%%
figure;
gscatter(score(:,1), score(:,2), labels_for_UMAP, [], [], [], 'off');
xlabel('Principal Component 1');
ylabel('Principal Component 2');
title('PCA of Feature Matrix');

% Optionally, add a legend if the labels are categorical
if iscategorical(labels_for_UMAP)
    legend(categories(labels_for_UMAP), 'Location', 'best');
else
    legend(arrayfun(@num2str, unique(labels_for_UMAP), 'UniformOutput', false), 'Location', 'best');
end

%%
disp('Saving')
timenow = char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss'));
save(['/media/user/teamshare/nick/behavior/grooming/', timenow, '_umap_dFF.mat'], 'data_for_UMAP', 'labels_for_UMAP', '-v7.3')

disp('Done')