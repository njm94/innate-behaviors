clear, clc, close all
addpath(fix_path('C:\Users\user\Documents\Nick\grooming\utils'))
data_root = 'Y:\nick\behavior\grooming\1p';
mice = {'ECR2_thy1', 'GER2_ai94', 'HYL3_tTA', 'IBL2_tTA'};

for j = 1:length(mice)
    load(fix_path([data_root, filesep, mice{j}, filesep, 'mask.mat']))
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


%%


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

all_maps = catcell(3, binary_maps);
% all_maps = catcell(3, all_behavior_maps);
% newnanmask = 1-isnan(mean(all_maps,3));
% newnanmask(newnanmask==0) = nan;
dicemat = zeros(size(all_maps,3), size(all_maps,3));
for i = 1:size(all_maps,3)
    for j = 1:size(all_maps,3)
        dicemat(j,i) = dice(all_maps(:,:,j), all_maps(:,:,i));
        % im1 = all_maps(:,:,j).*newnanmask;
        % im1 = (im1-min(im1(:)))./(max(im1(:))-min(im1(:)));
        % im1(isnan(im1)) = 0;
        % im2 = all_maps(:,:,i).*newnanmask;
        % im2 = (im2-min(im2(:)))./(max(im2(:))-min(im2(:)));
        % im2(isnan(im2)) = 0;
        % dicemat(j,i) = ssim(im1, im2, 'exponents', [1 1 1]);
    end
end

figure, 
imagesc(dicemat), 
% caxis([0 1]), 
% colormap gray
colormap(flipud(colormap(gray)))
colorbar


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



%%
% for i = 1:length(vars)
%     ax = figure(i);
%     box on
%     exportgraphics(ax, fix_path(['Y:\nick\behavior\grooming\figures\', char(vars(i)), '_dff_contours.png']), 'Resolution', 300)
% 
% end


%% compute pairwise dice
for i = 1:length(vars)
    trialcount = 1;
    for j = 1:length(mice)
        for k = 1:length(mice)
            if k > j
                dsim(trialcount, i) = dice(a{i}(:,:,j), a{i}(:,:,k));
                jsim(trialcount, i) = jaccard(a{i}(:,:,j), a{i}(:,:,k));
                trialcount = trialcount + 1;
            end
        end
    end
end




[p,tbl,stats] = anova1(dsim);
[c,m,h,gnames] = multcompare(stats);

figure, boxplot(dsim, 'Colors', 'k'),
xticklabels(["Lick", "Right", "Left", "Elliptical", "Right Asymmetric", ...
    "Left Asymmetric", "Bilateral"])
ylabel('Pairwise Dice Similarity Coefficient')
ax = gca;
ax.FontSize = 12;

title(['One-Way ANOVA, p = ', num2str(p, 2)])


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



%% single trial dice matrix

clear, clc


fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = '';
trial = [];
mousecount = 0;

left = nan(128, 128, length(data_list{1}));
right = nan(128, 128, length(data_list{1}));
elliptical = nan(128, 128, length(data_list{1}));
largeleft = nan(128, 128, length(data_list{1}));
largeright = nan(128, 128, length(data_list{1}));
largebilateral = nan(128, 128, length(data_list{1}));
lick = nan(128, 128, length(data_list{1}));

for j = 1:length(data_list{1})
    data_dir = data_list{1}{j};
    disp(['Starting ' data_dir])

    if isunix % working on linux computer - modify paths
        data_dir = strrep(data_dir, '\', '/');
        data_dir = strrep(data_dir, 'Y:/', '/media/user/teamshare/');
    end    
    [mouse_root_dir, exp_date, ~] = fileparts(data_dir);
    [~, mouse_id, ~] = fileparts(mouse_root_dir);
    [expt_root_dir, ~, ~] = fileparts(mouse_root_dir);

    new_mouse = ~strcmp(mouse_id, current_mouse);
    if new_mouse       
        load([mouse_root_dir filesep 'mask.mat'])
        mask = double(mask);
        mask(mask==0) = nan;
        trialcount = 1;
        mousecount = [mousecount mousecount(end)+1];
    else
        trialcount = trialcount + 1;
        mousecount = [mousecount mousecount(end)];
    end
    trial = [trial, trialcount];
    
    dff_path = [data_dir filesep 'outputs'];
    h = openfig([dff_path, filesep, getAllFiles(dff_path, '_dFF.fig')]);
    current_mouse = mouse_id;

    vars = {'left', 'right', 'lick', 'elliptical', 'largeleft', 'largeright', 'largebilateral'};

    for i = 1:length(vars)
        for jj = 1:length(h.Children)
            if strcmp(h.Children(jj).Title.String, vars{i})
                switch vars{i}
                    case 'left'
                        left(:,:,j) = h.Children(jj).Children.CData.*mask;
                    case 'right'
                        right(:,:,j) = h.Children(jj).Children.CData.*mask;
                    case 'lick'
                        lick(:,:,j) = h.Children(jj).Children.CData.*mask;
                    case 'elliptical'
                        elliptical(:,:,j) = h.Children(jj).Children.CData.*mask;
                    case 'largeleft'
                        largeleft(:,:,j) = h.Children(jj).Children.CData.*mask;
                    case 'largeright'
                        largeright(:,:,j) = h.Children(jj).Children.CData.*mask;
                    case 'largebilateral'
                        largebilateral(:,:,j) = h.Children(jj).Children.CData.*mask;
                end
            end
        end
    end


    close(h)

end
mousecount = mousecount(2:end);

%%


combined_bstack = cat(3, lick, left, right, elliptical, largeleft, largeright, largebilateral);
for i = 1:size(combined_bstack, 3)
    tmp = combined_bstack(:,:,i);
    thresh = prctile(tmp(:), 90);
    combined_bstack(:,:,i) = tmp > thresh;
end

%%

for i = 1:size(combined_bstack,3)
    if sum(combined_bstack(:,:,i), [1 2]) == 0, continue; end
    for j = 1:size(combined_bstack,3)   
        if sum(combined_bstack(:,:,j), [1 2]) == 0, continue; end
        dice_matrix(i,j) = dice(combined_bstack(:,:,i), combined_bstack(:,:,j));
    end
end


%%



figure, imagesc(dice_matrix)
colormap(bluewhitered)
yticks(25:25:175)
xticks(25:25:175)


%%


data_root = 'Y:\nick\behavior\grooming\1p';
mice = {'ECR2_thy1', 'GER2_ai94', 'HYL3_tTA', 'IBL2_tTA'};

for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    dff_path = [data_root, filesep, mice{j}, filesep, 'outputs'];
    h = openfig([dff_path, filesep, getAllFiles(dff_path, '_dFF.fig')]);
    rightmove(:,:,j) = h.Children(8).Children.CData;
    leftmove(:,:,j) = h.Children(10).Children.CData;
    left(:,:,j) = h.Children(16).Children.CData;
    right(:,:,j) = h.Children(14).Children.CData;
    lick(:,:,j) = h.Children(12).Children.CData;
    elliptical(:,:,j) = h.Children(24).Children.CData;
    largeleft(:,:,j) = h.Children(22).Children.CData;
    largeright(:,:,j) = h.Children(20).Children.CData;
    dropright(:,:,j) = h.Children(2).Children.CData;
    bilateral(:,:,j) = h.Children(18).Children.CData;
    close(h)
end

%% overlay contours from diff mice
clc
nanmask = zeros(128, 128, length(mice));

for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    nanmask(:,:,j) = mask;
end
nanmask(nanmask==0) = nan;





%%
clc, clear a v
tmp = [];

thresh = 80;
figure, axis off, hold on
for j = 1:length(mice)

    vars = ["lick", "right", "left", "elliptical", "largeright", "largeleft", "bilateral"];
         
    for i = 1:length(vars)    
        test = eval(vars(i));
        test = test .* nanmask;
        test = test(:,:,j);

        v = prctile(test(:), thresh);
        tmp = [tmp v];
        disp(v)
        a{i}(:,:,j) = test >= v;

        subplot(1,length(vars), i), axis off, hold on
            if v > 0
                contourf(flipud(test), [v v], 'FaceAlpha', 0.25)
    
                title(vars(i));
            end
%         end
    end
%     legend(labs, 'Location', 'Best')
end

%% compute pairwise dice
for i = 1:length(vars)
    trialcount = 1;
    for j = 1:length(mice)
        for k = 1:length(mice)
            if k > j
                similarity(trialcount, i) = dice(a{i}(:,:,j), a{i}(:,:,k));
                trialcount = trialcount + 1;
            end
        end
    end
end




[p,tbl,stats] = anova1(similarity);
[c,m,h,gnames] = multcompare(stats);

figure, boxplot(similarity, 'Colors', 'k'),
xticklabels(vars)
ylabel('Pairwise Dice Similarity Coefficient')
ax = gca;
ax.FontSize = 12;

title(['One-Way ANOVA, p = ', num2str(p, 2)])




%% not using


%%   do the same thing on ridge residuals

clear, clc


fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
data_list = textscan(fileID, formatSpec);

    left_counter = 1;
    right_counter = 1;
    elliptical_counter = 1;
    lick_counter = 1;
    largeleft_counter = 1;
    largeright_counter = 1;
    largebilateral_counter = 1;


for j = 1:length(data_list{1})%1:length(data_list{1})+1
    % load experiment specific data into cell array
    ridge_dir = [data_list{1}{j}, filesep, 'ridge_outputs_video'];
    if ~isfolder(ridge_dir)
        disp('Ridge on Videp not performed yet for current experiment date')
        continue
    end
    
    h = openfig([ridge_dir, filesep, getAllFiles(ridge_dir, 'residuals')]);

    vars = {'left', 'right', 'lick', 'elliptical', 'largeleft', 'largeright', 'largebilateral'};

    for i = 1:length(vars)
        for jj = 1:length(h.Children)
            if strcmp(h.Children(jj).Title.String, vars{i})
                switch vars{i}
                    case 'left'
                        left(:,:,left_counter) = h.Children(jj).Children.CData;
                        left_counter = left_counter + 1;
                    case 'right'
                        right(:,:,right_counter) = h.Children(jj).Children.CData;
                        right_counter = right_counter + 1;
                    case 'lick'
                        lick(:,:,lick_counter) = h.Children(jj).Children.CData;
                        lick_counter = lick_counter + 1;
                    case 'elliptical'
                        elliptical(:,:,elliptical_counter) = h.Children(jj).Children.CData;
                        elliptical_counter = elliptical_counter + 1;
                    case 'largeleft'
                        largeleft(:,:,largeleft_counter) = h.Children(jj).Children.CData;
                        largeleft_counter = largeleft_counter + 1;
                    case 'largeright'
                        largeright(:,:,largeright_counter) = h.Children(jj).Children.CData;
                        largeright_counter = largeright_counter + 1;
                    case 'largebilateral'
                        largebilateral(:,:,largebilateral_counter) = h.Children(jj).Children.CData;
                        largebilateral_counter = largebilateral_counter + 1;
                end
            end
        end
    end

    close(h)

end


%%

thresh = 95;
figure, hold on
for j = 1:length(vars)
    subplot(1, length(vars), j), hold on
    tmp0 = eval(vars{j});
%     tmp = mean(tmp0, 3);
    v = prctile(tmp(:), thresh);
%     imagesc(flipud(tmp)), colorbar

%     contourf(flipud(tmp), [v v], 'FaceAlpha', 0.25)
    for i = 1:size(tmp0,3)
%         tmp = abs(tmp0(:,:,i));
        tmp = -tmp0(:,:,i);
        v = prctile(tmp(:), thresh);
        contourf(flipud(tmp), [v v], 'FaceAlpha', 0.25)
        
    end
    title(vars{j})
end

%%
clear, clc

data_root = 'Y:\nick\behavior\grooming\1p';
mice = {'ECR2_thy1', 'GER2_ai94', 'HYL3_tTA', 'IBL2_tTA'};



for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    dff_path = [data_root, filesep, mice{j}, filesep, 'outputs'];
    h = openfig([dff_path, filesep, getAllFiles(dff_path, 'ridge_summary.fig')]);
    rightmove(:,:,j) = h.Children(2).Children.CData;
    leftmove(:,:,j) = h.Children(4).Children.CData;
    left(:,:,j) = h.Children(10).Children.CData;
    right(:,:,j) = h.Children(8).Children.CData;
    lick(:,:,j) = h.Children(6).Children.CData;
    elliptical(:,:,j) = h.Children(18).Children.CData;
    largeleft(:,:,j) = h.Children(16).Children.CData;
    largeright(:,:,j) = h.Children(14).Children.CData;
    bilateral(:,:,j) = h.Children(12).Children.CData;
    close(h)
end

%% overlay contours from diff mice
clc
nanmask = zeros(128, 128, length(mice));

for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    nanmask(:,:,j) = mask;
end
nanmask(nanmask==0) = nan;





%%
clc, clear a v
tmp = [];

thresh = 80;
figure, axis off, hold on
for j = 1:length(mice)

    vars = ["lick", "right", "left", "elliptical", "largeright", "largeleft", "bilateral"];
         
    for i = 1:length(vars)    
        test = eval(vars(i));
        test = test .* nanmask;
        test = test(:,:,j);

        v = prctile(test(:), thresh);
        a{i}(:,:,j) = test >= v;

        subplot(1,length(vars), i), axis off, hold on
            if v > 0
                contourf(flipud(test), [v v], 'FaceAlpha', 0.25)
    
                title(vars(i));
            end
%         end
    end
%     legend(labs, 'Location', 'Best')
end

%% compute pairwise dice
for i = 1:length(vars)
    trialcount = 1;
    for j = 1:length(mice)
        for k = 1:length(mice)
            if k > j
                similarity(trialcount, i) = dice(a{i}(:,:,j), a{i}(:,:,k));
                trialcount = trialcount + 1;
            end
        end
    end
end




[p,tbl,stats] = anova1(similarity);
[c,m,h,gnames] = multcompare(stats);

figure, boxplot(similarity, 'Colors', 'k'),
xticklabels(vars)
ylabel('Pairwise Dice Similarity Coefficient')
ax = gca;
ax.FontSize = 12;

title(['One-Way ANOVA, p = ', num2str(p, 2)])






%%   do the same thing on ridge unique explained var

clear, clc

data_root = 'Y:\nick\behavior\grooming\1p';
mice = {'ECR2_thy1', 'GER2_ai94', 'HYL3_tTA', 'IBL2_tTA'};



for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    dff_path = [data_root, filesep, mice{j}, filesep, 'outputs'];
    h = openfig([dff_path, filesep, getAllFiles(dff_path, 'ridge_summary.fig')]);
    rightmove(:,:,j) = h.Children(2).Children.CData;
    leftmove(:,:,j) = h.Children(4).Children.CData;
    left(:,:,j) = h.Children(10).Children.CData;
    right(:,:,j) = h.Children(8).Children.CData;
    lick(:,:,j) = h.Children(6).Children.CData;
    elliptical(:,:,j) = h.Children(18).Children.CData;
    largeleft(:,:,j) = h.Children(16).Children.CData;
    largeright(:,:,j) = h.Children(14).Children.CData;
    bilateral(:,:,j) = h.Children(12).Children.CData;
    close(h)
end

%% overlay contours from diff mice
clc
nanmask = zeros(128, 128, length(mice));

for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    nanmask(:,:,j) = mask;
end
nanmask(nanmask==0) = nan;





%%
clc, clear a v
tmp = [];

thresh = 80;
figure, axis off, hold on
for j = 1:length(mice)

    vars = ["lick", "right", "left", "elliptical", "largeright", "largeleft", "bilateral"];
         
    for i = 1:length(vars)    
        test = eval(vars(i));
        test = test .* nanmask;
        test = test(:,:,j);

        v = prctile(test(:), thresh);
        a{i}(:,:,j) = test >= v;

        subplot(1,length(vars), i), axis off, hold on
            if v > 0
                contourf(flipud(test), [v v], 'FaceAlpha', 0.25)
    
                title(vars(i));
            end
%         end
    end
%     legend(labs, 'Location', 'Best')
end

%% compute pairwise dice
for i = 1:length(vars)
    trialcount = 1;
    for j = 1:length(mice)
        for k = 1:length(mice)
            if k > j
                similarity(trialcount, i) = dice(a{i}(:,:,j), a{i}(:,:,k));
                trialcount = trialcount + 1;
            end
        end
    end
end




[p,tbl,stats] = anova1(similarity);
[c,m,h,gnames] = multcompare(stats);

figure, boxplot(similarity, 'Colors', 'k'),
xticklabels(vars)
ylabel('Pairwise Dice Similarity Coefficient')
ax = gca;
ax.FontSize = 12;

title(['One-Way ANOVA, p = ', num2str(p, 2)])

