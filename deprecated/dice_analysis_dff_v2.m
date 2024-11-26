clc, clear

addpath('/home/user/Documents/grooming/utils')
addpath('C:\Users\user\Documents\Nick\grooming\utils')

fileID = fopen('expt1_datalist.txt','r');

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
data_list{1} = fix_path(data_list{1});
current_mouse = '';

fs = 90 ;

thy1_idx = 1:7;
ai94_idx = 8:13;
hyl3_idx = 14:19;
ibl2_idx = 20:25;

spon_index = false(length(data_list{1}),1);
spon_index([1,2,8,9,14,15,20,21]) = true;

save_average_across_days = true;

load(fix_path('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat'));
% 
% bvars = ["DropLeft", "DropRight", "DropCenter", "Audio", ...
%         "Lick", ...
%         "Right", "Left", "Right Asymmetric", "Left Asymmetric" ...
%         "Elliptical", "Elliptical Left", "Elliptical Right" ...
%         "LeftMove", "RightMove"];


%%

for j = 1:length(data_list{1})
    % ignore spontaneous trials
    if spon_index(j), continue; end


    mouse_root = fileparts(data_list{1}{j});
    load([mouse_root, filesep, 'mask.mat'])
    load([mouse_root, filesep, 'atlas_tform.mat'])

    dff_path = [data_list{1}{j}, filesep, 'outputs'];
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

%%

right(right==0) = nan;
left(left==0) = nan;
ellip_left(ellip_left==0) = nan;
ellip_right(ellip_right==0) = nan;
largeright(largeright==0) = nan;
largeleft(largeleft ==0) = nan;
elliptical(elliptical==0)=nan;
lick(lick==0)=nan;


%%
clear avg_signal binary_map labs
% close all
vars = ["right", "left", "ellip_left", "ellip_right", "largeright", "largeleft", "elliptical", "lick"];
thresh = 90;
labcount = 1;
figure
for i = 1:length(vars)
    behavior_var = eval(vars(i));
    bcount = 1;
    for j = 1:size(behavior_var,3)
        tmp = behavior_var(:,:,j);
        % if no behaviors in the session, continue
        if isnan(mean(tmp, [1 2], 'omitnan')), continue; end

        avg_signal{i}(bcount, 1) = mean(tmp, [1 2], 'omitnan');

        v = prctile(tmp(:), thresh);
        binary_map{i}(:,:,bcount) = tmp >= v;
        data_map{i}(:,:,bcount) = tmp;
        labs{labcount} = vars(i);
        labcount = labcount + 1;
        bcount = bcount+1;
    end
    subplot(2,round(length(vars)/2), i),
    imagesc(mean(binary_map{i},3))
    caxis([0 1]), colormap(bluewhitered())

    axis off, hold on, set(gca, 'YDir', 'reverse');
    title(vars(i));
    for p = 1:length(dorsalMaps.edgeOutline)
        plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k', 'LineWidth', 1);
        xticks([])
        yticks([])
    end
end



%%
ax = gcf;
saveas(ax, fix_path(['Y:\nick\behavior\grooming\figures\','dFF_90p', '.svg']))
%%

avg_consolidated = uboxplot(cat(1, avg_signal{1}, avg_signal{2}), ...
    cat(1,avg_signal{3}, avg_signal{4}), ...
    cat(1, avg_signal{5}, avg_signal{6}), ...
    avg_signal{7}, avg_signal{8});


clc
figure
avg_dff = uboxplot(avg_signal{1}, avg_signal{2}, avg_signal{3}, avg_signal{4}, avg_signal{5}, avg_signal{6}, avg_signal{7}, avg_signal{8});
close()
[p,tbl,stats] = anova1(avg_dff);
[c,m,h,gnames] = multcompare(stats);

figure, boxplot(avg_dff, 'Colors', 'k'),
xticklabels(["Unilateral", "Elliptical Asymmetric", "Asymmetric", "Elliptical", "Lick"])
ylabel('Average \DeltaF/F_0')
ax = gca;
ax.FontSize = 12;

title(['One-Way ANOVA, p = ', num2str(p, 2)])
%%
% Create a table with the data, labeling each behavior
tbl = array2table(avg_consolidated, 'VariableNames', {'Unilateral', 'Elliptical Asymmetric', 'Asymmetric','Elliptical', 'Lick'}); 
% tbl = array2table(avg_dff, 'VariableNames', cellstr(vars)); 

% Define the repeated measures model
rm = fitrm(tbl, 'Unilateral-Lick ~ 1', 'WithinDesign', table([1:size(avg_consolidated,2)]', 'VariableNames', {'Behavior'}));
% rm = fitrm(tbl, [char(vars(1)),'-',char(vars(end)),' ~ 1'], 'WithinDesign', table([1:size(avg_dff,2)]', 'VariableNames', {'Behavior'}));

% Run repeated measures ANOVA
ranovaResults = ranova(rm);

% Display the results
disp(ranovaResults);

% If you want to conduct post hoc tests, you can use the multcompare function
% for pairwise comparisons between behaviors:
stats = multcompare(rm, 'Behavior');
%% compute pairwise dice
test = catcell(3, binary_map);
% test = catcell(3, data_map);
% newnanmask = 1-isnan(mean(test,3));
% newnanmask(newnanmask==0) = nan;
dicemat = zeros(size(test,3), size(test,3));

for i = 1:size(test,3)
%     disp(i)
    for j = 1:size(test,3)
        % im1 = test(:,:,j).*newnanmask;
        % im1(isnan(im1)) = nanmean(im1, [1 2]);
        % im1 = (im1-min(im1(:)))./(max(im1(:))-min(im1(:)));
        % im2 = test(:,:,i).*newnanmask;
        % im2(isnan(im2)) = nanmean(im2, [1 2]);
        % im2 = (im2-min(im2(:)))./(max(im2(:))-min(im2(:)));
        % dicemat(j,i) = ssim(im1, im2);
        dicemat(j,i) = dice(test(:,:,j), test(:,:,i));

        
    end
end
%%

figure, imagesc(dicemat), colormap(flipud(colormap(gray)))

%% pull out individual mice and look at similarity

vars = ["lick", "right", "left", "elliptical", "largeright", "largeleft", ...
    "ellip_right", "ellip_left"];
ecr2 = [];

for i = 1:length(vars)
    tmp = eval(vars(i));
    ecr2 = cat(3, ecr2, tmp(:,:,1:7));
end

nanidx = isnan(squeeze(mean(ecr2, [1 2], 'omitnan')));
num_behaviors = 
ecr2(:,:,nanidx) = [];

%%
dicemat = zeros(size(ecr2,3),size(ecr2,3));
for i =1 :size(ecr2,3)
    for j = 1:size(ecr2,3)
        im1 = ecr2(:,:,i);
        im2 = ecr2(:,:,j);

        v = prctile(im1(:), 80);
        im1 = im1>v;

        v = prctile(im2(:), 80);
        im2 = im2>v;
        dicemat(j,i) = dice(im1, im2);
    end
end

%%

figure, 
imagesc(dicemat), hold on
for i = 1:length(vars)-1
    hline(5*(i)+0.5, 'r-')
    vline(5*(i)+0.5, 'r-')
end
xticks(2.5:5:size(dicemat,1))
yticks(2.5:5:size(dicemat,1))
xticklabels(vars)
yticklabels(vars)
% caxis([0 1]), 
% colormap gray
colormap(flipud(colormap(gray)))
c = colorbar;
c.Label.String = 'Dice Similarity Coefficient';
title('DFF binary')



%% deprecated

% get the pairwise dice values from within each behavior (squares along the
% diagonal)
tmp_dicemat = tril(dicemat, -1);
tmp_dicemat(tmp_dicemat == 0) = nan;

count = 1;
dice_by_behavior = cell(size(vars));
for i = 1:length(binary_map)
    subset_dicemat = tmp_dicemat(count:count+size(binary_map{i},3)-1,count:count+size(binary_map{i},3)-1);
    disp(size(subset_dicemat))
    dice_by_behavior{i} = subset_dicemat(:);
    count = count + size(binary_map{i},3);
end

dsim = uboxplot(dice_by_behavior{1}, dice_by_behavior{2}, ...
    dice_by_behavior{3},dice_by_behavior{4}, dice_by_behavior{5}, ...
    dice_by_behavior{6}, dice_by_behavior{7}, dice_by_behavior{8});

[p,tbl,stats] = anova1(dsim);
[c,m,h,gnames] = multcompare(stats);

figure, boxplot(dsim, 'Colors', 'k'),
xticklabels(vars)
ylabel('Pairwise Dice Similarity Coefficient')
ax = gca;
ax.FontSize = 12;

title(['One-Way ANOVA, p = ', num2str(p, 2)])

%%
figure,
imagesc(imadjust(dicemat)), colormap(flipud(colormap(gray)));
hold on
counter = 0.5;
for i = 1:length(binary_map)
    counter = counter + size(binary_map{i},3);
    
    if i ~= length(binary_map)
        hline(counter, 'r-')
        vline(counter, 'r-')
    end
    new_x(i) = counter - size(binary_map{i},3)/2;
    new_y(i) = counter - size(binary_map{i},3)/2;
end
c = colorbar;
c.Label.String = 'Dice Similarity Coefficient'
xticks(new_x)
yticks(new_y)
xticklabels(vars)
yticklabels(vars)

%%

%%

test = catcell(3, binary_map);
% test(isnan(test)) = 0;
test = reshape(test, size(test,1)*size(test,2), [])';
Z = linkage(test, 'average', 'euclidean');


figure
[H, T, outperm] = dendrogram(Z);  % Plot the dendrogram
% [H, T, outperm] = dendrogram(Z,'Labels', labs, 'ColorThreshold', cutoff);  % Plot the dendrogram
% ylabel('Distance (1-correlation)')
xticks([])

%%
exportgraphics()


%%
% Create a table with the data, labeling each behavior
tbl = array2table(dsim, 'VariableNames', {'Lick', 'Right', 'Left', 'Elliptical', 'RightAsymmetric', 'LeftAsymmetric', 'EllipticalRight', 'EllipticalLeft'}); 

% Define the repeated measures model
rm = fitrm(tbl, 'Lick-EllipticalLeft ~ 1', 'WithinDesign', table([1:size(dsim,2)]', 'VariableNames', {'Behavior'}));

% Run repeated measures ANOVA
ranovaResults = ranova(rm);

% Display the results
disp(ranovaResults);

% If you want to conduct post hoc tests, you can use the multcompare function
% for pairwise comparisons between behaviors:
stats = multcompare(rm, 'Behavior');



%%   do the same thing on ridge unique explained var

clear, clc

data_root = fix_path('Y:\nick\behavior\grooming\1p');
mice = {'ECR2_thy1', 'GER2_ai94', 'HYL3_tTA', 'IBL2_tTA'};



for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
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
clear all_behavior_maps, close all
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

        subplot(2,round(length(vars)/2), i), axis off, hold on
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
% all_maps = catcell(3, binary_maps);
all_maps = catcell(3, all_behavior_maps);
newnanmask = 1-isnan(mean(all_maps,3));
newnanmask(newnanmask==0) = nan;
dicemat = zeros(size(all_maps,3), size(all_maps,3));
for i = 1:size(all_maps,3)
    for j = 1:size(all_maps,3)
%         dicemat(j,i) = dice(all_maps(:,:,j), all_maps(:,:,i));
        im1 = all_maps(:,:,j).*newnanmask;
        im1 = (im1-min(im1(:)))./(max(im1(:))-min(im1(:)));
        im1(isnan(im1)) = 0;
        im2 = all_maps(:,:,i).*newnanmask;
        im2 = (im2-min(im2(:)))./(max(im2(:))-min(im2(:)));
        im2(isnan(im2)) = 0;
        dicemat(j,i) = ssim(im1, im2, 'exponents', [1 1 1]);
    end
end

figure, 
imagesc(dicemat), 
% caxis([0 1]), 
% colormap gray
colormap(flipud(colormap(gray)))
colorbar

%%


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

