% Analyze signal decline over long periods of grooming 
% This is for 1p only


%% average events

clc, clear


if ~isunix
    addpath('C:\Users\user\Documents\Nick\ridgeModel');
    addpath('C:\Users\user\Documents\Nick\ridgeModel\widefield')
    addpath('C:\Users\user\Documents\Nick\ridgeModel\smallStuff') 
    addpath('C:\Users\user\Documents\Nick\grooming\utils')
    load('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat');

else
    addpath('/home/user/Documents/grooming/ridgeModel');
    addpath('/home/user/Documents/grooming/ridgeModel\widefield')
    addpath('/home/user/Documents/grooming/ridgeModel\smallStuff') 
    addpath('/home/user/Documents/grooming/utils')
    load('allenDorsalMap.mat');
end




fileID = fopen('expt1_datalist.txt','r');

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
cluster_data = fix_path('/media/user/teamshare/TM_Lab/nick/behavior/grooming/20241114092737_behavior_clustering.mat');
include_boris = true;
current_mouse = '';

fs = 90 ; % sampling rate
[b, a] = butter(2, 0.01/(fs/2), 'high');
aggregation_sz = 3; % window of time to aggregate behaviors
clen = 10; % duration after which continuous behavior is considered long
blen = 5; % duration of baseline period before and after continuous grooming event

thy1_idx = 1:7;
ai94_idx = 8:13;
camk_idx = 14:25;

save_average_across_days = true;

%%
 
for j = 4%1:length(data_list{1})+1
    disp('\n\n')
     try

        data_dir = fix_path(data_list{1}{j});
        disp(['Starting ' data_dir])

        [mouse_root_dir, exp_date, ~] = fileparts(data_dir);
        [~, mouse_id, ~] = fileparts(mouse_root_dir);
        [expt_root_dir, ~, ~] = fileparts(mouse_root_dir);
    
        new_mouse = ~strcmp(mouse_id, current_mouse);
    catch
        new_mouse = true;
    end
    if new_mouse && j<=length(data_list{1})     
        count = 1;
        disp('Loading master basis set')
        load([mouse_root_dir filesep 'Umaster.mat'])
        load([mouse_root_dir filesep 'mask.mat'])
        nanmask = mask;
        nanmask(nanmask==0) = nan;
        atlas_tform = load([mouse_root_dir filesep 'atlas_tform.mat']);
        clear Vmaster left right elliptical large_left large_right largebilateral lick LeftMove RightMove Audio drop_left drop_right
        current_mouse = mouse_id;

        fPath = [mouse_root_dir filesep 'outputs' filesep];

        behavior_frames = cell(1, 17);
    elseif new_mouse && j>length(data_list{1})
        disp('Finished')
        break
    end

    try
        brain_file = [data_dir, filesep, get_file_with_str(data_dir, 'cam0_svd')];
        beh_file = [data_dir, filesep, get_file_with_str(data_dir, [exp_date, '_svd'])];
        dlc_speed_file = [data_dir, filesep, get_file_with_str(data_dir, 'speed.csv')];
        timestamp_file = [data_dir, filesep, get_file_with_str(data_dir, 'trim.txt')];
        % snippets_dir = [data_dir, filesep, 'snippets'];
        boris_file = [data_dir, filesep, get_file_with_str(data_dir, 'events.tsv')];
    catch
        disp('Missing one or more of the required files. Skipping...')
        continue
    end

    %% load grooming events from BORIS file
    [events, b_idx, ~, trial_length, cluster_labels] = get_labels_from_clustering_results(cluster_data, boris_file, include_boris);

    
    % remove point events from event matrix
    % idx = contains(events.Properties.VariableNames, 'Drop') | ...
    %     contains(events.Properties.VariableNames, 'Video') | ...
    %     contains(events.Properties.VariableNames, 'Flail') ;
    idx = contains(events.Properties.VariableNames, 'Lick') | ...
        contains(events.Properties.VariableNames, 'Drop') | ...
        contains(events.Properties.VariableNames, 'Video') | ...
        contains(events.Properties.VariableNames, 'Flail') ;
    groom_events = removevars(events, idx);
    labels = groom_events.Properties.VariableNames;

    idx = contains(events.Properties.VariableNames, 'Video') | ...
        contains(events.Properties.VariableNames, 'Flail') ;
    eth_events = removevars(events, idx);
    
    bmat = any(table2array(groom_events),2);
    if ~any(bmat)
        continue
    end

    [episodes, idx] = aggregate(bmat, aggregation_sz);
    episode_durations = diff(idx, 1, 2)/fs;

    if ~any(episode_durations >= clen)
%         disp(['No episodes longer than ', num2str(clen),'s. Skipping...'])
        continue
    else
        idx = idx(episode_durations >= clen, :);
        if size(idx,1) > 1      
            [~, aa] = sort(diff(idx, 1, 2), 'descend');
            idx = idx(aa,:);
        end
        
    end

    % load experiment specific data into cell array
    load([data_dir filesep 'tform.mat'])


    disp('Loading brain data...')
    load(brain_file);
    U = permute(reshape(U, 128, 128, []), [2 1 3]);
    U = evaluate_tform(U, tformEstimate); % apply registration
    U = reshape(U, 128*128, []);

    % light-OFF signal may be captured in the brain data, resulting in a
    % massive filter artifact - remove a few frames just in case
    V = V(:,1:size(V, 2) - 5);

    % Brain video may have been resampled to wrong number of frames due to
    % inaccuracy in pre-processing step. Fix this by checking the length of
    % brain data and comparing to manually labeled VideoEnd event from
    % BORIS file. If there is a discrepancy, resample the brain data to the
    % length of the BORIS file to fix the issue.
    if size(V,2) ~= trial_length
        V = resamplee(V', trial_length, size(V,2))';
    end

    Vbrain = recastV(Umaster, U, s, V(:, 1:trial_length));
    Vmaster = filtfilt(b, a, Vbrain')';
   
    disp('Computing DF/F0')
%     dataR = permute(reshape(Umaster*Vmaster, 128, 128, []), [2 1 3]);
    dataR = reshape(Umaster*Vmaster, 128, 128, []);
    dFF = nanmask.*zscore((dataR-min(dataR(:))./mean(dataR - min(dataR(:)), 3))-1, [], 3);
    clear dataR Vbrain U s V

    %% load DLC tracks
    disp('Loading DLC tracks')
    dlc_speed = readmatrix(dlc_speed_file);
    fll_speed = dlc_speed(1:trial_length,1);
    flr_speed = dlc_speed(1:trial_length,2);
    
    % binarize movement speed
    fll_speed = fll_speed > mean(fll_speed) + std(fll_speed);
    flr_speed = flr_speed >  mean(flr_speed) + std(flr_speed);
    LeftMove = fll_speed;
    RightMove = flr_speed;
    
    clear dlc_speed fll_speed flr_speed

    % get ROIs - seeds are in atlas coords
    [seeds, labels] = get_seeds();

    any_behavior = any([bmat, LeftMove, RightMove], 2);
    tmp = any_behavior;

    for i = 1:size(idx,1)
        % crop data frames to include rest frames before and after grooming
        % episodes

        tmp_idx(i,:) = [idx(i,1)-round(blen*fs), idx(i,2)+round(blen*fs)];
        try 

            dFF_crop = dFF(:,:,tmp_idx(i,1):tmp_idx(i,2));

            % set windows
            prewin = 1:round(blen*fs);
            earlywin = round(blen*fs)+1:2*round(blen*fs);
            latewin = size(dFF_crop,3)-2*round(blen*fs):size(dFF_crop,3)-round(blen*fs);
            postwin = size(dFF_crop,3)-round(blen*fs)+1:size(dFF_crop,3);

            % fname = ['dFF_', mouse_id, '_', exp_date, '_', num2str(i)];
            % fileID = fopen(['/home/user/teamshare/TM_Lab/nick/behavior/grooming/1p/outputs/', filesep, fname, '.raw'],'w', 'l');
            % fwrite(fileID, single(dFF_crop),'single');
            % fclose(fileID);
    
            % get into atlas coords
            dFF_crop = imwarp(dFF_crop, atlas_tform.tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)), 'FillValues', nan);
            ts = getTimeseries(dFF_crop, seeds, 2);
            
            % if i == 1
            %     figure
            %     mask_warp = imwarp(mask, atlas_tform.tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)), 'FillValues', nan);
            %     imagesc(mask_warp), hold on
            %     scatter(seeds(:,1), seeds(:,2))
            % end
            % 
            % pre_corrmat{j}(:,:,i) = corrcoef(ts(prewin, :));
            % early_corrmat{j}(:,:,i) = corrcoef(ts(earlywin, :));
            % late_corrmat{j}(:,:,i) = corrcoef(ts(latewin, :));
            % post_corrmat{j}(:,:,i) = corrcoef(ts(postwin, :));
            % 
            % % continue
            % 
            % global_signal{j}{i,:} = squeeze(mean(dFF_crop, [1 2], 'omitnan'));
            % roi_signal{j}{i} = ts;
    
    
            % find a period of time that is equal in length with no activity
            groom_dur{j}(i) = size(ts,1)*fs - 2*round(blen*fs);

            % get behavior ethogram data
            eth_states{j}{i} = eth_events(tmp_idx(i,1):tmp_idx(i,2), :);

        catch
            continue
        end

    end



end


%% 
% states = ["Start", "Right", "Left", "Elliptical", ...
%     "Right Asymmetric", "Left Asymmetric", ...
%     "Elliptical Right", "Elliptical Left", "Stop", "Lick", "Drop"];
% close all
% ex_idx = [16 1; 5 1; 17 4; 19 2];
% for i = 1:size(ex_idx,1)
%     plot_ethogram(eth_states{ex_idx(i,1)}{ex_idx(i,2)}, states, fs);
%     axis([5 60, ylim])
% end







%%
% clc
% 7 2
% 4 1
didx = 7;
sidx = 2;
figure('Position', [186 304 1399 562])
states = ["Start", "Right", "Left", "Elliptical", ...
    "Right Asymmetric", "Left Asymmetric", ...
    "Elliptical Right", "Elliptical Left", "Stop", "Lick", "Drop"];
plot_ethogram(eth_states{didx}{sidx}, states, fs)
offset = 3.5;
t = xt(global_signal{didx}{sidx},fs,1);
hold on, plot(xt(global_signal{didx}{sidx},fs,1), global_signal{didx}{sidx}-offset, 'k')

% calculate slope of line excluding baseline period before and after
t_slope = t(round(blen*fs):end-round(blen*fs));
ex_slope = polyfit(t_slope, global_signal{didx}{sidx}(round(blen*fs):end-round(blen*fs)), 1);
plot(t_slope, t_slope.*ex_slope(1)-offset+ex_slope(2), 'r-')
line([0 t(end)], [-offset -offset], 'Color', [0 0 0], 'LineWidth', 1)
% line([1 1], [-1.5 -0.5], 'Color', [0 0 0], 'LineWidth', 2)
axis tight
yticks(-5.5:8.5)
yticklabels('auto')

% saveas(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'groom_episode_dff.svg']))




%%

clear mat_global sort_idx all_global all_dur

counter = 1;
all_dur = [];
all_groom = {};
for i = 1:length(global_signal)
    if ~isempty(global_signal{i})
        for j = 1:length(global_signal{i})
            if ~isempty(global_signal{i}{j})
                % all_global{counter} = global_signal{i}{j}(1:end-(round(blen*fs)));
                all_global{counter} = global_signal{i}{j};
                all_dur(counter) = length(all_global{counter});
                counter = counter + 1;
            end
        end
    end
end

maxdur = max(all_dur);
[~,sort_idx] = sort(all_dur);
mat_global = nan(length(all_global), maxdur);
for i = 1:length(all_global)
    mat_global(i, 1:length(all_global{sort_idx(i)})) = all_global{sort_idx(i)};
end

% hold on
% vline(0, 'r-')

%%

% append a row and column of zeros since pcolor ignores the last row and column
figure, hold on
t = xt([mat_global(1,:), 0], fs, 2)-5;
num_ep = size(mat_global,1)+1:-1:1;
pcolor(t, num_ep, padarray(mat_global, [1, 1], 0, 'post')), 
colormap(bluewhitered())
shading flat
caxis([-2 8])
line([40 50], [50 50], 'Color', [0 0 0], 'LineWidth', 2)
xticklabels([])
yticklabels([])
% set(gca, 'YDir', 'reverse')
axis([t(1) t(end) 1 size(mat_global,1)+1])
line([0 0], [0 size(mat_global,1)+1], 'Color', [0 0 0], 'LineWidth', 2)
% colorbar
p50 = prctile(all_dur, 50);
line([t(p50) t(p50)], [0 size(mat_global,1)+1], 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '--')
box on
xticks([])
yticks([])

% exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'longgroom_all.png']), 'Resolution', 300)
%%
figure, axis off, colorbar, caxis([-2 8]), colormap(bluewhitered())
exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'longgroom_colorbar.png']), 'Resolution', 300)
%%

figure, hold on
shadedErrorBar(t(1:p50), mean(mat_global(:,1:p50), 'omitnan'), std(mat_global(:,1:p50), 'omitnan'), 'k', 1)
line([0 0], ylim, 'Color', [0 0 0], 'LineWidth', 2)

xticklabels([])
yticklabels([])
exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'longgroom_average.png']), 'Resolution', 300)

%%
t = xt(mat_global, fs, 2)-5;
t_pre = t<=0;
% t_on = abs(t)<=1;
t_early = t>0 & t<=5;
dff_pre = mean(mat_global(:, t_pre), 2, 'omitnan');
% dff_on = mean(mat_global(:, t_on), 2, 'omitnan');
dff_early = mean(mat_global(:,t_early), 2, 'omitnan');


for i = 1:length(all_dur)
    t_late_idx = all_dur(sort_idx(i));
    t_late = t_late_idx - 2*blen*fs;
    dff_late(i) = mean(mat_global(i, t_late:t_late_idx-blen*fs), 'omitnan');
    dff_post(i) = mean(mat_global(i, t_late_idx-blen*fs + 1: end), 'omitnan');
    % t_late = t_late_idx - blen*fs;
    % dff_late(i) = mean(mat_global(i, t_late:end), 'omitnan');
end
% t_late = t>5;
% dff_late = mean(mat_global(:,t_late), 2, 'omitnan');


plot_data = [dff_pre, dff_early, dff_late', dff_post'];
figure, boxplot(plot_data, 'Colors', 'k', 'Symbol', '')
hold on
swarmchart(repmat([1 2 3 4], size(dff_early,1), 1), plot_data, 'k', 'XJitterWidth', 0.25)
ylabel('Mean cortical \DeltaF/F_0 (\sigma)')
xticklabels({'Before 5s', 'First 5s', 'Last 5s', 'After 5s'})
ax = gca;
ax.FontSize = 14;


% [p,~,stats] = anova1([dff_on, dff_early, dff_late])
% [c,m,h,gnames] = multcompare(stats);

clc
% Create a table with the data
tbl = array2table(plot_data, 'VariableNames', {'Before', 'Early', 'Late', 'After'});

% Define the repeated measures model
rm = fitrm(tbl, 'Before-After ~ 1', 'WithinDesign', [1 2 3 4]'); % 3 time periods

% Run repeated-measures ANOVA
ranova_results = ranova(rm);

% Display the results
disp(ranova_results);

% Perform pairwise comparisons
pairwise_results = multcompare(rm, 'Time', 'ComparisonType', 'bonferroni');
disp(pairwise_results);

% saveas(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'long_groom_binnedstates.svg']))

%%
clear p
for i = 1:length(all_global)
    
    X = all_global{i}(round(blen*fs):end);
    t = xt(X, fs, 1);
    avg_slope(i,:) = polyfit(t, X, 1);
    
    % figure, plot(t, X)
    % f = fit(t,X,'exp1')
    % hold on, plot(t, f)
    % hfdhd
    % linslope(i) = p(1);
    % disp(p)
end
close all
figure, boxplot(avg_slope(:,1), 'Colors', 'k', 'Symbol', '')
hold on, hline(0, 'k--')
swarmchart(ones(size(avg_slope,1),1), avg_slope(:,1), 'k', 'XJitterWidth', 0.5)
ylabel('Slope of best-fit line')
xticks([])
ax = gca; 
ax.FontSize = 14;


[h, p, ci, stats] = ttest(avg_slope(:,1), 0, 'Tail', 'left');

% Display results
disp(['p-value: ', num2str(p)]);
disp(['Test statistic (t): ', num2str(stats.tstat)]);
disp(['95% confidence interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);

if h == 1
    disp('The mean is significantly less than zero.');
else
    disp('The mean is not significantly less than zero.');
end

% saveas(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'long_groom_slope.svg']))
%%
save('/media/user/teamshare/nick/behavior/grooming/1p/outputs/long_groom_global_5spre.mat', 'all_global')
%%
load(fix_path('/media/user/teamshare/nick/behavior/grooming/1p/ECR2_thy1/mask.mat'))
load(fix_path('/media/user/teamshare/nick/behavior/grooming/1p/ECR2_thy1/atlas_tform.mat'))
h = openfig(fix_path('/media/user/teamshare/nick/behavior/grooming/1p/ECR2_thy1/atlas_aligned.fig'));

mask = imwarp(mask, tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
mask(mask==0) = nan;

for i = 1:length(h.Children.Children)
    if strcmp(h.Children.Children(i).Type, 'line')
        h.Children.Children(i).LineWidth = 2;
        h.Children.Children(i).Color = [1 1 1];
    else
        h.Children.Children(i).CData = h.Children.Children(i).CData .* uint16(mask);
        set(h.Children.Children(i), 'AlphaData', ~isnan(mask))
    end
end

scalebar_length = 2000/dorsalMaps.desiredPixelSize;
line([400 400+scalebar_length], [520 520], 'Color', [0 0 0], 'LineWidth', 2)

% rois = {'MOS_1-L', 'SSP-ul-L', 'SSP-bfd-L', 'MOP_1-L', 'VIS-p-L', 'RSP_2-L'};
% hold on
% for i = 1:length(rois)
%     tmp = seeds(strcmp(labels, rois(i)),:);
%     plot(tmp(1), tmp(2), '.', 'MarkerSize', 30)
% end

%%
exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'example1pimage.png']), 'Resolution', 300)




%%





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
% figure('Position', [82 474 1811 128])
example_mask(example_mask==0)=nan;
for i = 1:length(vars)
    figure
    tmp = eval(vars(i));
    
    % subplot(1,length(vars), i)
    pcolor(example_mask.*tmp(:,:,example_mouse));
    shading interp
    set(gca, 'YDir', 'reverse');
    axis off, box off
    % imagesc(example_mask.*tmp(:,:,example_mouse));
    xticks([])
    yticks([])
    caxis([0 2.5])
    % title(vars(i))
    % colorbar
    % saveas(gcf, fix_path(['Y:\nick\behavior\grooming\figures\',char(vars(i)), '_example.pdf']))
    % exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\',char(vars(i)), '_example.png']), 'Resolution', 300)
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
close all
clear all_behavior_maps
thresh = 80;
vars = ["lick", "right", "left", "elliptical", "largeright", "largeleft", "ellip_right", "ellip_left"];
% figure, axis off, hold on
for j = 1:length(mice)
    load(fix_path([data_root, filesep, mice{j}, filesep, 'atlas_tform.mat']));
    for i = 1:length(vars)
        
        behavior_map = eval(vars(i));
        
        behavior_map = behavior_map(:,:,j);
        behavior_map = nanmask(:,:,j).*imwarp(behavior_map, tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
        all_behavior_maps{i}(:,:,j) = behavior_map;
        v = prctile(behavior_map(:), thresh);
        binary_maps{i}(:,:,j) = behavior_map >= v;

        % subplot(1,round(length(vars)), i), axis off, hold on
        figure(i), axis off, hold on
        for p = 1:length(dorsalMaps.edgeOutline)
            plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k', 'LineWidth', 1);
            xticks([])
            yticks([])
        end
            if v > 0
                contourf(binary_maps{i}(:,:,j).*j, [j-0.1 j-0.1], 'FaceAlpha', 0.25)

                % title(vars(i));
            else
                disp('v is not greater than 0')
                disp(vars(i))
            end
            set(gca, 'YDir', 'reverse');
    end
    
end

