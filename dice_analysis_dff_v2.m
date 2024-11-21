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



clc
figure
avg_dff = uboxplot(avg_signal{1}, avg_signal{2}, avg_signal{3}, avg_signal{4}, avg_signal{5}, avg_signal{6}, avg_signal{7}, avg_signal{8});
close()
[p,tbl,stats] = anova1(avg_dff);
[c,m,h,gnames] = multcompare(stats);

figure, boxplot(avg_dff, 'Colors', 'k'),
xticklabels(vars)
ylabel('Average \DeltaF/F_0')
ax = gca;
ax.FontSize = 12;

title(['One-Way ANOVA, p = ', num2str(p, 2)])
%%

%% compute pairwise dice
test = catcell(3, binary_map);
dicemat = zeros(size(test,3), size(test,3));

for i = 1:size(test,3)
%     disp(i)
    for j = 1:size(test,3)
        dicemat(j,i) = dice(test(:,:,j), test(:,:,i));
%         disp(j)
%         if isnan(test(i)) || isnan(test(j)), continue; end
        
    end
end

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
        hline(counter, 'r:')
        vline(counter, 'r:')
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
xticklabels(vars)
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

