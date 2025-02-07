% Analyze changes in functional connectivity associated with long duration
% continuous grooming behaviors
%
% Expected outcome - Increase in FC at onset, followed by reduced FC at the
% end
%
% This is on 1p data only




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
 
for j = 1:length(data_list{1})+1
    disp('\n\n')
     try

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
    catch
        new_mouse = true;
    end
    if new_mouse 
        if ~isempty(current_mouse)
            % if all the trials have completed for the previous mouse,
            % compile all behavior frames and take the average
%%
            if save_average_across_days
                figure
                for iii = 1:size(behavior_frames,2)
                    % transpose the mean image to get the brains in proper
                    % orientation
                    if isempty(behavior_frames{iii}), continue; end
                    mean_image = mean(behavior_frames{iii},3)';
                    subplot(4,4,iii)
                    imagesc(mean_image)
%                     set(gca, 'clim', [-max(abs(mean_image(:))) max(abs(mean_image(:)))])
                    c=colorbar;
                    c.Label.String = '\DeltaF/F_0 (\sigma)';
                    title(bvars(iii))
                    xticks([])
                    yticks([])
%                     colormap(bluewhitered())
%                     colormap(fliplr(redblue([],[],'w')))
                end
                                       
                savefig(gcf, [fPath,  char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_dFF.fig'])
            end
            
            %%
            % all mice completed - break the loop
            if j == length(data_list{1})+1, break; end
        end        
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

        
    end

    try
        brain_file = [data_dir, filesep, get_file_with_str(data_dir, 'cam0_svd')];
        beh_file = [data_dir, filesep, get_file_with_str(data_dir, [exp_date, '_svd'])];
%         ME_file = [data_dir, filesep, get_file_with_str(data_dir, 'MEsvd')];
%         frame_file = [data_dir, filesep, get_file_with_str(data_dir, 'singleFrame')];
%         grooming_file = [data_dir, filesep, 'grooming_events_roi_filtered.mat'];
        dlc_speed_file = [data_dir, filesep, get_file_with_str(data_dir, 'speed.csv')];
        timestamp_file = [data_dir, filesep, get_file_with_str(data_dir, 'trim.txt')];
        % snippets_dir = [data_dir, filesep, 'snippets'];
        boris_file = [data_dir, filesep, get_file_with_str(data_dir, 'events.tsv')];
    catch
        disp('Missing one or more of the required files. Skipping...')
        continue
    end

    %% load grooming events from BORIS file
    events = read_boris(boris_file);

    % create the aggregated behavior episodes - if none exceed 15s, skip
    % consolidate lick events
    lick_idx = contains(events.Properties.VariableNames, 'Lick');
    lick_events = events(:,lick_idx);
    lick_events = any(table2array(lick_events),2);
    
    % remove lick and point events from event matrix
    idx = contains(events.Properties.VariableNames, 'Lick') | ...
        contains(events.Properties.VariableNames, 'Drop') | ...
        contains(events.Properties.VariableNames, 'Video') | ...
        contains(events.Properties.VariableNames, 'Flail') ;
    stroke_events = removevars(events, idx);
    labels = stroke_events.Properties.VariableNames;

    bmat = any(table2array(stroke_events),2);
    [episodes, idx] = aggregate(bmat, aggregation_sz);
    episode_durations = diff(idx, 1, 2)/90;

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
    trial_length = size(V, 2) - 5; 
    V = V(:,1:trial_length);

    % Brain video may have been resampled to wrong number of frames due to
    % inaccuracy in pre-processing step. Fix this by checking the length of
    % brain data and comparing to manually labeled VideoEnd event from
    % BORIS file. If there is a discrepancy, resample the brain data to the
    % length of the BORIS file to fix the issue.
    trial_length = find(events.("Video End"));
    if size(V,2) ~= trial_length
        V = resamplee(V', trial_length, size(V,2))';
    end

    Vbrain = recastV(Umaster, U, s, V(:, 1:trial_length));
    Vmaster = filtfilt(b, a, Vbrain')';

    % transform one last time to get into coordinates of atlas
%     U = permute(reshape(U, 128, 128, []), [2 1 3]);
%     U = evaluate_tform(U, atlas_tform.tform); % apply registration
%     U = reshape(U, 128*128, []);
    
    disp('Computing DF/F0')
%     dataR = permute(reshape(Umaster*Vmaster, 128, 128, []), [2 1 3]);
    dataR = reshape(Umaster*Vmaster, 128, 128, []);
    dFF = zscore((dataR-min(dataR(:))./mean(dataR - min(dataR(:)), 3))-1, [], 3);
    clear dataR Vbrain U s V

    %% load DLC tracks
    disp('Loading DLC tracks')
    dlc_speed = readmatrix(dlc_speed_file);
    fll_speed = dlc_speed(1:trial_length,1);
    flr_speed = dlc_speed(1:trial_length,2);
    
    % binarize movement speed
    fll_speed = fll_speed > mean(fll_speed) + std(fll_speed);
    flr_speed = flr_speed >  mean(flr_speed) + std(flr_speed);
    
    k = 90; % smoothing kernel - tp ensure rest phase doesn't capture movement
    LeftMove= movmax(fll_speed, k);
    RightMove= movmax(flr_speed, k);

    

    clear dlc_speed fll_speed flr_speed

    % get ROIs - seeds are in atlas coords
    [seeds, labels] = get_seeds();

    % use these vars to find contiguous stretch of stationary behavior
    % which is same duration as grooming- we will overwrite tmp in the loop
    % starting with longest grooming behavior first, so we avoid choosing
    % same time
    any_behavior = any([bmat, lick_events, LeftMove, RightMove], 2);
    tmp = any_behavior;

    for i = 1:size(idx,1)
        % crop data frames to include rest frames before and after grooming
        % episodes

        tmp_idx(i,:) = [idx(i,1)-round(blen*fs), idx(i,2)+round(blen*fs)];
        try dFF_crop = dFF(:,:,tmp_idx(i,1):tmp_idx(i,2));
            

        
%         clear pc90
%         stepsize = round(0.25 * fs);
%         corrk = 1 * fs;
%         count = 1;
%         for ii = 1:stepsize:size(ts,1)-corrk
%             pc90(count) = pca_video_90(dFF_crop(:,:,ii:ii+corrk));
%             count = count +1;
%         end

    

        %%

        % get into atlas coords
        dFF_crop = imwarp(dFF_crop, atlas_tform.tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
        ts = getTimeseries(dFF_crop, seeds, 2);

        % set windows
        prewin = 1:round(blen*fs);
        earlywin = round(blen*fs)+1:2*round(blen*fs);
        latewin = size(ts,1)-2*round(blen*fs):size(ts,1)-round(blen*fs);
        postwin = size(ts,1)-round(blen*fs)+1:size(ts,1);

        pre_corrmat{j}(:,:,i) = corrcoef(ts(prewin, :));
        early_corrmat{j}(:,:,i) = corrcoef(ts(earlywin, :));
        late_corrmat{j}(:,:,i) = corrcoef(ts(latewin, :));
        post_corrmat{j}(:,:,i) = corrcoef(ts(postwin, :));

        % calculate complexity (num PCs to reach 99%) for each window
        dFF_crop = dFF_crop.* imwarp(nanmask, atlas_tform.tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)), 'FillValues', nan);
        dFF_crop = reshape(dFF_crop, size(dFF_crop,1)*size(dFF_crop,2), []);
        dFF_crop(isnan(dFF_crop(:,1)),:) = [];

        global_signal{j}{i} = mean(dFF_crop);
        roi_signal{j}{i} = ts;
        
        
        
        % [~, ~, latent] = pca(dFF_crop(:, prewin)', 'economy', true); 
        % % Compute the cumulative explained variance
        % explained_variance = cumsum(latent) / sum(latent);    
        % % Find the number of components that explain 99% of the variance
        % pre_numComponents{j}(i) = find(explained_variance >= 0.99, 1);
        % 
        % [~, ~, latent] = pca(dFF_crop(:, earlywin)', 'economy', true); 
        % explained_variance = cumsum(latent) / sum(latent);    
        % early_numComponents{j}(i) = find(explained_variance >= 0.99, 1);
        % 
        % [~, ~, latent] = pca(dFF_crop(:, latewin)', 'economy', true); 
        % explained_variance = cumsum(latent) / sum(latent);    
        % late_numComponents{j}(i) = find(explained_variance >= 0.99, 1);
        % 
        % [~, ~, latent] = pca(dFF_crop(:, postwin)', 'economy', true); 
        % explained_variance = cumsum(latent) / sum(latent);    
        % post_numComponents{j}(i) = find(explained_variance >= 0.99, 1);


%         % perform a moving window correlation over the data
        % stepsize = round(0.25 * fs);
        % corrk = 1 * fs;
        % count = 1;
        % for ii = 1:stepsize:size(ts,1)-corrk
        %     corrmat{j}{i}(:,:,count) = corrcoef(ts(ii:ii+corrk,:));
        %     count = count +1;
        % end

        % find a period of time that is equal in length with no activity
        groom_dur{j}(i) = size(ts,1)*fs - 2*round(blen*fs);
        % 
        % a = strfind(tmp', zeros(1,groom_dur));
        % if isempty(a)
        %     continue
        % else
        %     dFF_crop = dFF(:,:,a(1):a(1)+groom_dur);
        %     dFF_crop = imwarp(dFF_crop, atlas_tform.tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
        % 
        %     ts_rest = getTimeseries(dFF_crop, seeds, 2);
        %     rest_corrmat{j}(:,:,i) = corrcoef(ts_rest);
        % 
        %     tmp(a(1):a(1)+groom_dur)=1;
        %     count = count + 1;
        % end
        catch
            continue
        end

    end



end



%% starting here
midx = {'thy1_idx', 'ai94_idx', 'hyl3_idx', 'ibl2_idx'};
m = 4;


test_pre = pre_corrmat;
test_pre = test_pre(:,eval(midx{m}));
test_pre = catcell(3, test_pre(~cellfun(@isempty, test_pre)));



upper_t = triu(true(size(test_pre,1)), 1); % Mask for upper triangle
test_pre = reshape(test_pre(repmat(upper_t, 1, 1, size(test_pre,3))), [], size(test_pre,3));
test_pre = tanh(mean(atanh(test_pre)))';
% test_pre = mean(test_pre)';


test_early = early_corrmat;
test_early = test_early(:,eval(midx{m}));
test_early = catcell(3, test_early(~cellfun(@isempty, test_early)));
test_early = reshape(test_early(repmat(upper_t, 1, 1, size(test_early,3))), [], size(test_early,3));
test_early = tanh(mean(atanh(test_early)))';
% test_early = mean(test_early)';

test_late = late_corrmat;
test_late = test_late(:,eval(midx{m}));
test_late = catcell(3, test_late(~cellfun(@isempty, test_late)));
test_late = reshape(test_late(repmat(upper_t, 1, 1, size(test_late,3))), [], size(test_late,3));
test_late = tanh(mean(atanh(test_late)))';
% test_late = mean(test_late)';


test_post = post_corrmat;
test_post = test_post(:,eval(midx{m}));
test_post = catcell(3, test_post(~cellfun(@isempty, test_post)));
test_post = reshape(test_post(repmat(upper_t, 1, 1, size(test_post,3))), [], size(test_post,3));
test_post = tanh(mean(atanh(test_post)))';
% test_post = mean(test_post)';

% remove zeros for now but check where they are coming from in the future

if any(test_pre==0) 
    zidx=find(test_pre==0); 
    test_pre(zidx) = [];
    test_early(zidx) = [];
    test_late(zidx) = [];
    test_post(zidx) = [];
end

plot_data = [test_pre, test_early, test_late, test_post];

figure, boxplot(plot_data, 'Colors', 'k', 'Symbol', '')
hold on
swarmchart(repmat([1 2 3 4], size(test_early,1), 1), plot_data, 'k', 'XJitterWidth', 0.25)
ylabel('Average correlation')
xticklabels({'Before 5s', 'First 5s', 'Last 5s', 'After 5s'})
ax = gca;
ax.FontSize = 14;




% Run Friedman test
[p_friedman, tbl, stats] = friedman(plot_data, 1); % 1 indicates repeated measures
disp(['Friedman test p-value: ', num2str(p_friedman)]);

% If the Friedman test is significant, perform post-hoc pairwise comparisons
if p_friedman < 0.05
    fprintf('Post-hoc pairwise comparisons:\n');
    
    % Define time points
    time_labels = {'Pre', 'Early', 'Late', 'Post'};
    
    % Get all possible pairwise combinations
    pairs = nchoosek(1:4, 2);
    p_values = zeros(size(pairs, 1), 1);
    
    % Perform Wilcoxon signed-rank tests for each pair
    for i = 1:size(pairs, 1)
        [p_values(i), ~, stats] = signrank(plot_data(:, pairs(i, 1)), plot_data(:, pairs(i, 2)));
        fprintf('%s vs. %s: p = %.4f\n', time_labels{pairs(i, 1)}, time_labels{pairs(i, 2)}, p_values(i));
    end

    % Multiple comparisons correction (Bonferroni)
    corrected_p_values = min(p_values * length(p_values), 1); % Bonferroni correction
    fprintf('\nBonferroni corrected p-values:\n');
    
    for i = 1:size(pairs, 1)
        fprintf('%s vs. %s: p_corrected = %.4f\n', time_labels{pairs(i, 1)}, time_labels{pairs(i, 2)}, corrected_p_values(i));
    end
end

timepoints = {'test_pre', 'test_early', 'test_late', 'test_post'};
figure,  hold on
for i = 1:length(timepoints)
    [a,b] = ksdensity(eval(timepoints{i}));
    plot(a,b, 'LineWidth', 2);
end
legend({'Before 5s', 'First 5s', 'Last 5s', 'After 5s'}, 'Location', 'Best')



%%
hyl3_idx = 14:19;
ibl2_idx = 20:25;
midx = {'thy1_idx', 'ai94_idx', 'hyl3_idx', 'ibl2_idx'};
m = 3;
caxlims = [0.5 1];
test_pre = pre_corrmat;
test_pre = pre_corrmat(:,eval(midx{m}));
test_pre = catcell(3, test_pre(~cellfun(@isempty, test_pre)));
disp(size(test_pre,3))

test_early = early_corrmat;
test_early = early_corrmat(:,eval(midx{m}));
test_early = catcell(3, test_early(~cellfun(@isempty, test_early)));

test_late = late_corrmat;
test_late = late_corrmat(:,eval(midx{m}));
test_late = catcell(3, test_late(~cellfun(@isempty, test_late)));

test_post = post_corrmat;
test_post = post_corrmat(:,eval(midx{m}));
test_post = catcell(3, test_post(~cellfun(@isempty, test_post)));


figure
subplot(2,4,1), imagesc(tanh(nanmean(atanh(test_pre),3))), colorbar, caxis(caxlims)
xticks(1:size(test_pre,1))
yticks(1:size(test_pre,1))
xticklabels(labels)
yticklabels(labels)
subplot(2,4,2), imagesc(tanh(nanmean(atanh(test_early), 3))), colorbar, caxis(caxlims)
xticks(1:size(test_pre,1))
yticks(1:size(test_pre,1))
xticklabels(labels)
yticklabels(labels)
subplot(2,4,3), imagesc(tanh(nanmean(atanh(test_late), 3))), colorbar, caxis(caxlims)
xticks(1:size(test_pre,1))
yticks(1:size(test_pre,1))
xticklabels(labels)
yticklabels(labels)

subplot(2,4,4), imagesc(tanh(nanmean(atanh(test_post), 3))), colorbar, caxis(caxlims)
xticks(1:size(test_pre,1))
yticks(1:size(test_pre,1))
xticklabels(labels)
yticklabels(labels)

colormap(flipud(colormap(gray)))
%%

tb = tanh(mean(atanh(test_pre),3));
tb = tb(upper_t);

te = tanh(mean(atanh(test_early),3));
te = te(upper_t);

tl = tanh(mean(atanh(test_late),3));
tl = tl(upper_t);

ta = tanh(mean(atanh(test_post),3));
ta = ta(upper_t);


figure, plot(sort(tb), (1:numel(tb))/numel(tb), 'LineWidth', 2);
hold on
plot(sort(te), (1:numel(te))/numel(te), 'LineWidth', 2)
plot(sort(tl), (1:numel(tl))/numel(tl), 'LineWidth', 2)
plot(sort(ta), (1:numel(ta))/numel(ta), 'LineWidth', 2)
legend({'Before 5', 'First 5', 'Last 5', 'After 5'}, 'FontSize', 16, 'Location', 'Best')



%%
figure
d_early_pre = tanh(nanmean(atanh(test_early),3))-tanh(nanmean(atanh(test_pre), 3));
subplot(1,3,1), imagesc(d_early_pre), 
caxis([-0.2 0.2])
xticks(1:size(test_pre,1))
yticks(1:size(test_pre,1))
xticklabels(labels)
yticklabels(labels)
title('First 5s - Before 5s')
c=colorbar;
c.Label.String = '\Delta Correlation';

d_late_pre = tanh(nanmean(atanh(test_late),3))-tanh(nanmean(atanh(test_pre), 3));
subplot(1,3,2), imagesc(d_late_pre),  caxis([-0.2 0.2])
xticks(1:size(test_pre,1))
yticks(1:size(test_pre,1))
xticklabels(labels)
yticklabels(labels)
title('Last 5s - Before 5s')
c=colorbar;
c.Label.String = '\Delta Correlation';

d_post_pre = tanh(nanmean(atanh(test_post), 3))-tanh(nanmean(atanh(test_pre), 3));
subplot(1,3,3), imagesc(d_post_pre), colorbar, caxis([-0.2 0.2])
xticks(1:size(test_pre,1))
yticks(1:size(test_pre,1))
xticklabels(labels)
yticklabels(labels)
title('After 5s - Before 5s')
c=colorbar;
c.Label.String = '\Delta Correlation';

colormap(bluewhitered())
% for i = 1:length(data_list{1})
%     subplot(1,3,1)
%     imagesc()
% end



d_time_points = {'d_early_pre', 'd_late_pre', 'd_post_pre'};
time_points = {'test_early', 'test_late', 'test_post'};
% figure,

data_1 = atanh(test_pre);
% data_1 = test_pre;



for i = 1:length(time_points)
    d_tmp = eval(d_time_points{i});
    % data_2 = eval(time_points{i});
    data_2 = atanh(eval(time_points{i}));

    clear p_values
    for ii = 1:size(data_1,1)
        for jj = 1:size(data_2,1)
            [p_values(ii,jj), observeddifference, effectsize] = permutationTest(squeeze(data_1(ii,jj,:)), squeeze(data_2(ii,jj,:)), 1000); 
        end
    end

    % for ii = 1:size(data_1,1)
    %     for jj = 1:size(data_2,1)
    %         [~, p_values(ii,jj)] = ttest(squeeze(data_1(ii,jj,:)), squeeze(data_2(ii,jj,:))); 
    %     end
    % end


    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_values(upper_t),0.05,'dep','yes');
    h_tmp = zeros(size(p_values));
    h_tmp(upper_t) = h;

    if nansum(h(:)) ==0, h_tmp = zeros(size(d_tmp)); end
    G = graph(d_tmp.*h_tmp, 'upper');
    LWidths = abs(G.Edges.Weight)*20;
    LColors = zeros(numedges(G), 3);

    for ii = 1:numedges(G)
        if G.Edges.Weight(ii) <0
            LColors(ii,:) = [0 0 1];
        else
            LColors(ii,:) = [1 0 0];
        end
    end
    figure, hold on
    
    for p = 1:length(dorsalMaps.edgeOutline)
        plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k', 'LineWidth', 1);
        xticks([])
        yticks([])
    end
    set(gca, 'YDir', 'reverse');
    hold on, axis equal off
    h=plot(G, 'XData', seeds(:,1), 'YData', seeds(:,2), 'EdgeColor', LColors, 'LineWidth', LWidths, 'NodeLabel', {}, 'NodeColor', 'k');
end
%%







%%


iii = 5;
figure, imagesc(mean(test_post(:,:,25:end) - test_pre(:,:,25:end),3))
colorbar
colormap(bluewhitered())
%%
%%
grooming_bout = 3;
w = 15;
nanmask = mask;
nanmask(nanmask==0) = nan;
dFF_crop = dFF(:,:,idx(grooming_bout,1)-w*fs:idx(grooming_bout,2)+w*fs);
dFF_crop = dFF_crop .* nanmask;



t = xt(dFF_crop, fs, 3)-w;

figure, plot(t, squeeze(nanmean(dFF_crop, [1 2])), 'k', 'LineWidth', 1)
hold on

for i = 1:size(stroke_events,2)
    tmp = table2array(stroke_events(:,i));
    tmp = tmp(idx(grooming_bout,1)-w*fs:idx(grooming_bout,2)+w*fs);
    if any(tmp)
        if strcmpi(stroke_events.Properties.VariableNames{i}, 'Right') || strcmpi(stroke_events.Properties.VariableNames{i}, 'Left')
            col = 'm';
        else
            col = 'c';
        end
        patchplot(t(arr2idx(tmp)), [-2 3], col, 0.1)
    end
    % patchplot(t([w*fs size(dFF_crop,3)-(w*fs)]), [-2 3], 'c', 0.1)
end
axis tight

xlabel('Time (s)')
ylabel('\DeltaF/F_0 (\sigma)')

%%

%%

figure
pp = diff(behaviors{4},1,2);
for i = 1:size(behaviors{4},1)
    subplot(1,3,i)
    tmp = mean(dFF(:,:,behaviors{4}(i,1):behaviors{4}(i,2)),3);
    imagesc(tmp)
    title(num2str(pp(i)))
    yticks([])
    xticks([])
end

%%


counter = 1;
all_dur = [];
all_groom = {};
for i = 1:length(global_signal)
    if ~isempty(global_signal{i})
        for j = 1:length(global_signal{i})
            all_global{counter} = global_signal{i}{j}(1:end-(round(blen*fs)));
            all_dur(counter) = length(all_global{counter});
            counter = counter + 1;
        end
    end
end

maxdur = max(all_dur);
[~,sort_idx] = sort(all_dur);
mat_global = nan(length(all_global), maxdur);
for i = 1:length(all_global)
    mat_global(i, 1:length(all_global{sort_idx(i)})) = all_global{sort_idx(i)};
end

figure, 
t = xt(mat_global, fs, 2)-blen;
joyPlot(mat_global', t, 2, 'FaceColor', 'w', 'FaceAlpha', 1)
% hold on
% vline(0, 'r-')


%%

t_on = abs(t)<=1;
t_early = t>1 & t<=5;
t_late = t>5;

dff_on = mean(mat_global(:, t_on), 2, 'omitnan');
dff_early = mean(mat_global(:,t_early), 2, 'omitnan');
dff_late = mean(mat_global(:,t_late), 2, 'omitnan');

plot_data = [dff_on, dff_early, dff_late];
figure, boxplot(plot_data, 'Colors', 'k', 'Symbol', '')
hold on
swarmchart(repmat([1 2 3], size(dff_on,1), 1), plot_data, 'k', 'XJitterWidth', 0.25)
ylabel('Mean cortical \DeltaF/F_0 (\sigma)')
xticklabels({'Onset [-1, 1]', 'Early [1, 5]', 'Late [5, end]'})
ax = gca;
ax.FontSize = 14;


% [p,~,stats] = anova1([dff_on, dff_early, dff_late])
% [c,m,h,gnames] = multcompare(stats);


% Create a table with the data
tbl = array2table(plot_data, 'VariableNames', {'Onset', 'Early', 'Late'});

% Define the repeated measures model
rm = fitrm(tbl, 'Onset-Late ~ 1', 'WithinDesign', [1 2 3]'); % 3 time periods

% Run repeated-measures ANOVA
ranova_results = ranova(rm);

% Display the results
disp(ranova_results);

% Perform pairwise comparisons
pairwise_results = multcompare(rm, 'Time', 'ComparisonType', 'bonferroni');
disp(pairwise_results);

%%
clear p
for i = 1:length(all_global)
    
    X = all_global{i}(round(blen*fs):end);
    t = xt(X, fs, 2);
    p(i,:) = polyfit(t, X, 1);
    
    % figure, plot(t, X)
    % f = fit(t,X,'exp1')
    % hold on, plot(t, f)
    % hfdhd
    % linslope(i) = p(1);
    % disp(p)
end
close all
figure, boxplot(p(:,1), 'Colors', 'k', 'Symbol', '')
hold on, hline(0, 'k--')
swarmchart(ones(size(p,1),1), p(:,1), 'k', 'XJitterWidth', 0.5)
ylabel('Slope of best-fit line')
xticks([])
ax = gca; 
ax.FontSize = 14;


[h, p, ci, stats] = ttest(p(:,1), 0, 'Tail', 'left');

% Display results
disp(['p-value: ', num2str(p)]);
disp(['Test statistic (t): ', num2str(stats.tstat)]);
disp(['95% confidence interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);

if h == 1
    disp('The mean is significantly less than zero.');
else
    disp('The mean is not significantly less than zero.');
end
%%
save('/media/user/teamshare/nick/behavior/grooming/1p/outputs/long_groom_global_5spre.mat', 'all_global')

%%
figure
% test = randn(1000, 10);
joyPlot(ts, xt(ts,fs,1), 0.5, 'FaceColor', 'w', 'FaceAlpha', 1)

%%


% %%
%     figure('Position', [152 426 768 496])
%     for ii = 1:length(labels)
%         groom_data = [];
%         prewin = round(90*3);
%         postwin = round(90*3);
% 
%         tmp = eval(labels{ii});
%         if ~isempty(tmp)
%             for i = 1:size(tmp,1)
%                 groom_data(:,:,:,i) = dFF(:,:,tmp(i,1) - prewin:tmp(i,1) + postwin);
%             end
%             
%             subplot(length(grooming_type), 1, ii)
%             imagesc(mean([groom_data(:,:,1:90), groom_data(:,:,91:180), ...
%                 groom_data(:,:,181:270), groom_data(:,:,271:360), ...
%                 groom_data(:,:,361:450), groom_data(:,:,451:540)],3)), colorbar
%             colormap(colormap_blueblackred());
%             c=colorbar;
%             c.Label.String = '\DeltaF/F_0 (z-score)';
%             caxis([-4 4])
%             yticks([])
%             xticks(64:128:128*6)
%             xticklabels({'-2.5', '-1.5', '-0.5', '0.5', '1.5', '2.5'})
%             xlabel('Time from grooming onset (s)')
%             title(grooming_type{ii})
%         end
%     end
%     ax = gcf;
% %     exportgraphics(ax,[mouse_root_dir, filesep,exp_date,'.pdf'],'ContentType','vector')
% 
% 
% end