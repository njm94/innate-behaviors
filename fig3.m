%%
clear, clc, close all
addpath(fix_path('C:\Users\user\Documents\Nick\grooming\utils'))
data_root = 'Y:\nick\behavior\grooming\1p';
mice = {'ECR2_thy1', 'GER2_ai94', 'HYL3_tTA', 'IBL2_tTA'};
example_mouse = 1;
for j = 1:length(mice)
    load(fix_path([data_root, filesep, mice{j}, filesep, 'mask.mat']))
    if j == example_mouse
        example_mask = mask;
    end
    dff_path = fix_path([data_root, filesep, mice{j}, filesep, 'outputs']);
    dff_fig = getAllFiles(dff_path, '_dFF.fig');
    % If there are multiple versions, use the most recent one
    if size(dff_fig,1) > 1
        dff_fig = sort(dff_fig);
        dff_fig = dff_fig{end};
    end
    h = openfig([dff_path, filesep, dff_fig]);
    for i = 1:length(h.Children)
        switch(h.Children(i).Title.String)
            case 'RightMove' 
                rightmove(:,:,j) = h.Children(i).Children(end).CData;
            case 'LeftMove' 
                leftmove(:,:,j) = h.Children(i).Children(end).CData;
            case 'largebilateral'
                bilateral(:,:,j) = h.Children(i).Children(end).CData;
            case 'Elliptical'
                elliptical(:,:,j) = h.Children(i).Children(end).CData;
            case 'Elliptical Right'
                ellip_right(:,:,j) = h.Children(i).Children(end).CData;
            case 'Elliptical Left'
                ellip_left(:,:,j) = h.Children(i).Children(end).CData;
            case 'Left Asymmetric'
                largeleft(:,:,j) = h.Children(i).Children(end).CData;
            case 'Right Asymmetric'
                largeright(:,:,j) = h.Children(i).Children(end).CData;
            case 'Left'
                left(:,:,j) = h.Children(i).Children(end).CData;
            case 'Right'
                right(:,:,j) = h.Children(i).Children(end).CData;
            case 'Lick'
                lick(:,:,j) = h.Children(i).Children(end).CData;
            otherwise
                continue
        end
    end

    close(h)
end

%%  plot example mouse dff
vars = ["lick", "right", "left", "elliptical", "largeright", "largeleft", "ellip_right", "ellip_left"];
figure
for i = 1:length(vars)
    
    tmp = eval(vars(i));
    subplot(1,length(vars), i)
    imagesc(example_mask.*tmp(:,:,example_mouse));
    xticks([])
    yticks([])
%     caxis([0 2.5])
    title(vars(i))
%     colorbar
end

%% overlay contours from diff mice
clc
% nanmask = zeros(128, 128, length(mice));

load('allenDorsalMap.mat');
clear nanmask
for j = 1:length(mice)
    load(fix_path([data_root, filesep, mice{j}, filesep, 'mask.mat']))
    atlas_tform = load(fix_path([data_root, filesep, mice{j}, filesep, 'atlas_tform.mat']));
    warpmask = imwarp(mask, atlas_tform.tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
    nanmask(:,:,j) = warpmask;
    % nanmask(:,:,j) = mask;
end
nanmask(nanmask==0) = nan;



%%
clear all_behavior_maps
thresh = 80;
vars = ["lick", "right", "left", "elliptical", "largeright", "largeleft", "ellip_right", "ellip_left"];
figure, axis off, hold on
for j = 1:length(mice)
    load(fix_path([data_root, filesep, mice{j}, filesep, 'atlas_tform.mat']));
    for i = 1:length(vars)
        
        behavior_map = eval(vars(i));
        
        behavior_map = behavior_map(:,:,j);
        behavior_map = nanmask(:,:,j).*imwarp(behavior_map, tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
        all_behavior_maps{i}(:,:,j) = behavior_map;
        v = prctile(behavior_map(:), thresh);
        binary_maps{i}(:,:,j) = behavior_map >= v;

        subplot(1,round(length(vars)), i), axis off, hold on
        for p = 1:length(dorsalMaps.edgeOutline)
            plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k', 'LineWidth', 1);
            xticks([])
            yticks([])
        end
            if v > 0
                contourf(binary_maps{i}(:,:,j).*j, [j-0.1 j-0.1], 'FaceAlpha', 0.25)

                title(vars(i));
            else
                disp('v is not greater than 0')
                disp(vars(i))
            end
            set(gca, 'YDir', 'reverse');
    end
end

%%

all_maps = catcell(3, binary_maps);
dicemat = zeros(size(all_maps,3), size(all_maps,3));
for i = 1:size(all_maps,3)
    for j = 1:size(all_maps,3)
        dicemat(j,i) = dice(all_maps(:,:,j), all_maps(:,:,i));
    end
end

%%
figure, 
imagesc(dicemat), hold on
for i = 1:length(vars)-1
    hline(length(mice)*(i)+0.5, 'r-')
    vline(length(mice)*(i)+0.5, 'r-')
end
xticks(2.5:length(mice):size(dicemat,1))
yticks(2.5:length(mice):size(dicemat,1))
xticklabels(vars)
yticklabels(vars)
% caxis([0 1]), 
% colormap gray
colormap(flipud(colormap(gray)))
c = colorbar;
c.Label.String = 'Dice Similarity Coefficient';
title('DFF binary')






%% check the betas for an example mouse

load('Y:\nick\behavior\grooming\1p\ECR2_thy1\outputs\2024-11-21-15-15-38_cvFull.mat')
load('Y:\nick\behavior\grooming\1p\ECR2_thy1\mask.mat')
load('Y:\nick\behavior\grooming\1p\ECR2_thy1\Umaster.mat')

%%
fullLabels
visual = true;
vars = ["Lick", "Right", "Left", "Elliptical", "Right Asymmetric", ...
    "Left Asymmetric", "Elliptical Right", "Elliptical Left"];
cBetaRight = check_beta('Audio', fullLabels, fullIdx, Umaster, fullBeta{1}, Vfull, [], visual);

%%
figure
for i = 1:length(vars)
    cBetaRight = check_beta(vars(i), fullLabels, fullIdx, Umaster, fullBeta{1}, Vfull, [], visual);
    mbeta(:,:,i) = mask.*peak2peak(cBetaRight,3);
    subplot(1,length(vars), i)
    imagesc(mbeta(:,:,i)), 
    colorbar%caxis([0 30]),
    title(vars(i))
    xticks([])
    yticks([])
end

%%   do the same thing on ridge unique explained var

clear, clc

data_root = fix_path('Y:\nick\behavior\grooming\1p');
mice = {'ECR2_thy1', 'GER2_ai94', 'HYL3_tTA', 'IBL2_tTA'};
example_mouse = 1;

vars = ["Lick", "right", "left", "elliptical", "largeright", "largeleft", "ellip_right", "ellip_left"];
for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    if j == example_mouse
        example_mask = mask;
    end
    dff_path = [data_root, filesep, mice{j}, filesep, 'outputs'];
    ridge_fig = getAllFiles(dff_path, 'summary.fig');
    % If there are multiple versions, use the most recent one
    if size(ridge_fig,1) > 1
        ridge_fig = sort(ridge_fig);
        ridge_fig = ridge_fig{end};
    end
    disp(ridge_fig)

    h = openfig([dff_path, filesep, ridge_fig]);
    for i = 1:length(h.Children)
        switch(h.Children(i).Title.String)
            case 'RightMove' 
                rightmove(:,:,j) = h.Children(i).Children(end).CData;
            case 'LeftMove' 
                leftmove(:,:,j) = h.Children(i).Children(end).CData;
            case 'largebilateral'
                bilateral(:,:,j) = h.Children(i).Children(end).CData;
            case 'Elliptical'
                elliptical(:,:,j) = h.Children(i).Children(end).CData;
            case 'Elliptical Right'
                ellip_right(:,:,j) = h.Children(i).Children(end).CData;
            case 'Elliptical Left'
                ellip_left(:,:,j) = h.Children(i).Children(end).CData;
            case 'Left Asymmetric'
                largeleft(:,:,j) = h.Children(i).Children(end).CData;
            case 'Right Asymmetric'
                largeright(:,:,j) = h.Children(i).Children(end).CData;
            case 'Left'
                left(:,:,j) = h.Children(i).Children(end).CData;
            case 'Right'
                right(:,:,j) = h.Children(i).Children(end).CData;
            case 'Lick'
                lick(:,:,j) = h.Children(i).Children(end).CData;
            otherwise
                continue
        end
    end
    close(h)
end

%% overlay contours from diff mice
clc
% nanmask = zeros(128, 128, length(mice));

load('allenDorsalMap.mat');
clear nanmask
for j = 1:length(mice)
    load(fix_path([data_root, filesep, mice{j}, filesep, 'mask.mat']))
    atlas_tform = load(fix_path([data_root, filesep, mice{j}, filesep, 'atlas_tform.mat']));
    warpmask = imwarp(mask, atlas_tform.tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
    nanmask(:,:,j) = warpmask;
    % nanmask(:,:,j) = mask;
end
nanmask(nanmask==0) = nan;





%%
clear all_behavior_maps, %close all
thresh = 80;
vars = ["lick", "right", "left", "elliptical", "largeright", "largeleft", "ellip_right", "ellip_left"];
figure, axis off, hold on
for j = 1:length(mice)
    load(fix_path([data_root, filesep, mice{j}, filesep, 'atlas_tform.mat']));
    for i = 1:length(vars)
        
        behavior_map = eval(vars(i));
        
        behavior_map = behavior_map(:,:,j);
        behavior_map = nanmask(:,:,j).*imwarp(behavior_map, tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
        all_behavior_maps{i}(:,:,j) = behavior_map;
        v = prctile(behavior_map(:), thresh);
        binary_maps{i}(:,:,j) = behavior_map >= v;

        subplot(1,round(length(vars)), i), axis off, hold on
        for p = 1:length(dorsalMaps.edgeOutline)
            plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k', 'LineWidth', 1);
            xticks([])
            yticks([])
        end
            if v > 0
                contourf(binary_maps{i}(:,:,j).*j, [j-0.1 j-0.1], 'FaceAlpha', 0.25)

                title(vars(i));
            else
                disp('v is not greater than 0')
                disp(vars(i))
            end
            set(gca, 'YDir', 'reverse');
    end
end

%%
all_maps = catcell(3, binary_maps);
dicemat = zeros(size(all_maps,3), size(all_maps,3));
for i = 1:size(all_maps,3)
    for j = 1:size(all_maps,3)
        dicemat(j,i) = dice(all_maps(:,:,j), all_maps(:,:,i));
    end
end
%%
figure, 
imagesc(dicemat), hold on
for i = 1:length(vars)-1
    hline(length(mice)*(i)+0.5, 'r-')
    vline(length(mice)*(i)+0.5, 'r-')
end
xticks(2.5:length(mice):size(dicemat,1))
yticks(2.5:length(mice):size(dicemat,1))
xticklabels(vars)
yticklabels(vars)
% caxis([0 1]), 
% colormap gray
colormap(flipud(colormap(gray)))
c = colorbar;
c.Label.String = 'Dice Similarity Coefficient';
title('Unique explained variance')




%% compare dice similarity matrix within behavior vs between behavior
clear within_subset between_subset
% within behavior
N = length(mice);
for i = 1:length(vars)
    widx = (i-1)*N;
    tmp = dicemat(widx+(1:N),widx+(1:N));
    tmp(logical(triu(tmp))) = nan;
    within_subset(:,i) = tmp(~isnan(tmp));

    tmp = dicemat(widx+(1:N), :);
    tmp(:, widx+(1:N)) = nan;
    between_subset(:,i) = tmp(~isnan(tmp));
end

figure, 
for i = 1:length(vars)
    subplot(length(vars), 1, i)
    histogram(between_subset(:,i), 20), hold on
    vline(mean(within_subset(:,i)), 'k-')
    xlim([0 1])
    title(vars(i))
end


% figure, 
for i = 1:length(vars)
    % subplot(length(vars), 1, i)
    % histogram(between_subset(:,i), 20), hold on
    tmp = mean(within_subset(:,i));
    my_p(i) = mean(tmp>between_subset(:,i));
    % vline(mean(within_subset(:,i)), 'k-')
    % xlim([0 1])
    % title(vars(i))
end

