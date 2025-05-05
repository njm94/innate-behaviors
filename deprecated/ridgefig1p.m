clear, close all, clc


% cd('/media/user/teamshare/nick/behavior/grooming/code/')
addpath('/media/user/teamshare/nick/behavior/grooming/code/ridgeModel')
addpath('/media/user/teamshare/nick/behavior/grooming/code/ridgeModel/widefield')
addpath('/media/user/teamshare/nick/behavior/grooming/code/ridgeModel/smallStuff')

% addpath(fix_path('C:\Users\user\Documents\Nick\ridgeModel'))
% addpath(fix_path('C:\Users\user\Documents\Nick\ridgeModel\widefield'))
% addpath(fix_path('C:\Users\user\Documents\Nick\ridgeModel\smallStuff'))
% addpath(fix_path('C:\Users\user\Documents\Nick\grooming\utils'))

%%

fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
groom_list = textscan(fileID, formatSpec);

fileID = fopen('expt3_datalist.txt','r');
formatSpec = '%s';
lick_list = textscan(fileID, formatSpec);

% data_list = [groom_list{1}; lick_list{1}];
data_list = fix_path(groom_list{1});
current_mouse = '';


behaviors = {'Audio', 'Drop Left', 'Drop Right', 'Drop Center', 'Lick', ...
    'Elliptical', 'Elliptical Right', 'Elliptical Left', 'Left', ...
    'Left Asymmetric', 'Right', 'Right Asymmetric', 'FLL', 'FLR'};

grooming_behaviors = {'Elliptical', 'Elliptical Right', ...
    'Elliptical Left', 'Left', 'Left Asymmetric', 'Right', ...
    'Right Asymmetric'};


thy1_idx = 1:7;
ai94_idx = 8:13;
hyl3_idx = 14:19;
ibl2_idx = 20:25;

trial_to_analyze = 7;
if trial_to_analyze < ai94_idx(1)
    mouse = 'thy1';
    mcounter = 0;
elseif trial_to_analyze < hyl3_idx(1)
    mouse = 'ai94';
    mcounter = 7;
elseif trial_to_analyze < ibl2_idx(1)
    mouse = 'hyl3';
    mcounter = 13;
else
    mouse = 'ibl2';
    mcounter = 19;
end


% grooming labels from unsupervised clustering
cluster_data = fix_path('Y:\nick\behavior\grooming\20241114092737_behavior_clustering.mat');
include_boris = true;
fs = 90;


%%
data_dir = data_list{trial_to_analyze};
[mouse_root_dir, exp_date, ~] = fileparts(data_dir);
[~, mouse_id, ~] = fileparts(mouse_root_dir);
[expt_root_dir, ~, ~] = fileparts(mouse_root_dir);

disp('Loading master basis set, mask')
load([mouse_root_dir filesep 'Umaster.mat'])
load([mouse_root_dir filesep 'mask.mat'])

disp('Loading ridge variables')
cvfull_path = getAllFiles([mouse_root_dir, filesep, 'outputs'], 'cvFull.mat');
% if there are multiple files, choose the last one
if size(cvfull_path, 1) > 1, cvfull_path = cvfull_path{end}; end
load([mouse_root_dir, filesep, 'outputs', filesep, cvfull_path])

switch mouse
    case 'thy1'
        num_trials = 7;
    otherwise
        num_trials = 6;
end

trial_length = size(Vfull,2)/num_trials;
sub_idx = trial_length*(trial_to_analyze-1-mcounter)+1:trial_length*(trial_to_analyze-mcounter-1)+trial_length;

disp('Load experiment specific brain data')
brain_file = [data_dir, filesep, getAllFiles(data_dir, 'cam0_svd')];
boris_file = [data_dir, filesep, getAllFiles(data_dir, 'events.tsv')];
    
load([data_dir filesep 'tform.mat'])
load(brain_file);
U = permute(reshape(U, 128, 128, []), [2 1 3]);
U = evaluate_tform(U, tformEstimate); % apply registration
U = reshape(U, 128*128, []);

% recast V, filter, resample to behavior, then crop to length of Vmodel
Vbrain = recastV(Umaster, U, s, V(:, 1:size(V,2)-5));
[b, a] = butter(2, 0.01/(fs/2), 'high');
Vbrain = filtfilt(b, a, Vbrain')';


disp('Reading BORIS')
[events, b_idx, ~, video_end, cluster_labels] = get_labels_from_clustering_results(cluster_data, boris_file, include_boris);

% resample brain to match video_end size
Vbrain = resamplee(Vbrain', video_end, size(Vbrain,2))';
Vbrain = Vbrain(:, 1:trial_length);

%%

disp('Reconstructing model and real data')
data_modeled = reshape(Umaster * Vfull(:,sub_idx), [128, 128, trial_length]);
data_real = reshape(Umaster * Vbrain, [128, 128, trial_length]);

%% get seeds for evaluating model fit
clear u v
load([mouse_root_dir filesep 'atlas_tform.mat'])

template_img = loadtiff(fix_path([mouse_root_dir, filesep, 'template.tif']));
[seeds, labels] = get_seeds();
close()
rois = {'MOS_1-L', 'SSP-ul-L', 'SSP-bfd-L', 'MOP_1-L', 'VIS-p-L', 'RSP_2-L'};
figure(1), 
imagesc(template_img), colormap gray,
hold on
for i = 1:length(rois)
    tmp = seeds(strcmp(labels, rois(i)),:);
    % seeds are in atlas coords. Use inverse of atlas tform to get them
    % back in original space

    [u(i), v(i)] = transformPointsInverse(tform, tmp(1), tmp(2));
    plot(u(i),v(i), '.', 'MarkerSize', 30)
end

%%
nanmask = mask;
nanmask(mask==0) = nan;
t = xt(data_real, fs, 3);
dff = data_real-min(data_real(:));
dff = zscore((dff - mean(dff, 3))./mean(dff,3), [], 3);
%%
figure('Position', [428 453 817 256]), 
plot(t, squeeze(mean(nanmask.*dff, [1 2], 'omitnan')),' k')
hold on, axis tight
groom_events = any(table2array(events(:,1:end-1)),2);
patchplot(t(arr2idx(groom_events(1:size(dff,3)))), [7.5 8], 'k', 1)
vline(912, 'r-')
vline(960, 'r-')
line([31 90], [4 4], 'Color', [0 0 0], 'LineWidth', 2)
line([31 31], [4 6], 'Color', [0 0 0], 'LineWidth', 2)
axis off
% saveas(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'groom_wholesession_dff.svg']))
%%
% Aggregate groom events, but this time use the modeled data since each 
% event has a surrounding window.

groom_events = arr2idx(squeeze(any(data_modeled, [1 2])));
clear expVar episode_duration

for i = 1:size(groom_events,1)
    % tmp_model = reshape(nanmask.*(data_modeled(:,:,groom_events(i,1):groom_events(i,2))), 128*128, []);
    % tmp_real = reshape(nanmask.*(data_real(:,:,groom_events(i,1):groom_events(i,2))), 128*128, []);

    episode_duration(i) = (groom_events(i,2)-groom_events(i,1))/fs;
    expVar(:,:,i) = reshape(modelCorr(Vbrain(:, groom_events(i,1):groom_events(i,2)), ...
        Vfull(:,sub_idx(1)-1+groom_events(i,1):sub_idx(1)-1+groom_events(i,2)), Umaster.*nanmask(:)), 128, 128);
end

%%

figure, scatter(episode_duration-2.5, squeeze(mean(expVar, [1 2], 'omitnan')))



%%
close all
clear u v
load([mouse_root_dir filesep 'atlas_tform.mat'])

template_img = loadtiff(fix_path([mouse_root_dir, filesep, 'template.tif']));
[seeds, labels] = get_seeds();
close()
rois = {'MOP_1-L', 'SSP-ul-L', 'PTLp-L'};
figure(1), 
imagesc(template_img), colormap gray,
hold on
axis off
for i = 1:length(rois)
    tmp = seeds(strcmp(labels, rois(i)),:);
    % seeds are in atlas coords. Use inverse of atlas tform to get them
    % back in original space

    [u(i), v(i)] = transformPointsInverse(tform, tmp(1), tmp(2));
    plot(u(i),v(i), '.', 'MarkerSize', 60)
end
% exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'example_map_with_rois.png']), 'Resolution', 300)

%%

betas_to_plot = {'Lick', 'Right', 'Elliptical', 'Right Asymmetric', 'Elliptical Right'};
fullLabels
visual = false;
cBetaRight = check_beta('Elliptical Right', fullLabels, fullIdx, Umaster, mean(catcell(3,fullBeta),3), Vfull, [], visual);
% close all
close gcf

for i = 1:length(betas_to_plot)
    tmpBeta = check_beta(betas_to_plot{i}, fullLabels, fullIdx, Umaster, mean(catcell(3,fullBeta),3), Vfull, [], visual);
    % figure(2, 'Position', [675 69 570 797])
    % subplot(length(betas_to_plot), 1, i)    
    imagesc([mask.*tmpBeta(:,:,45), mask.*tmpBeta(:,:,90), mask.*tmpBeta(:,:,135)]);
    caxis([-20 60])
    % hold on
    % for j = 1:length(u)
    %     line([u(j)-2 u(j)-2], [v(j)-2 v(j)+2], 'Color', [1 1 1], 'LineWidth', 2)
    %     line([u(j)+2 u(j)+2], [v(j)-2 v(j)+2], 'Color', [1 1 1], 'LineWidth', 2)
    %     line([u(j)-2 u(j)+2], [v(j)-2 v(j)-2], 'Color', [1 1 1], 'LineWidth', 2)
    %     line([u(j)-2 u(j)+2], [v(j)+2 v(j)+2], 'Color', [1 1 1], 'LineWidth', 2)
    % end
    % colorbar
    colormap(bluewhitered())
    axis equal off
    exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', betas_to_plot{i}, '_beta.png']), 'Resolution', 300)
    close gcf
    
    tmprois = getTimeseries(tmpBeta, [u; v]', 2);
    t = xt(tmprois, fs, 1);
    % figure(3)
    % subplot(length(betas_to_plot), 1, i)
    plot(t-0.5, tmprois)
    hold on
    axis off
    vline(0, 'k-')
    vline(0.5, 'k-')
    vline(1, 'k-')
    % saveas(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', betas_to_plot{i}, '_timeseries.svg']))
    % close gcf
end




%%
rois_real = getTimeseries(data_real, [u; v]', 2);
rois_modeled = getTimeseries(data_modeled, [u; v]', 2);

%%
cols = [[0 0.4470 0.7410];
    [0.8500 0.3250 0.0980];
    [0.9290 0.6940 0.1250];
    [0.4940 0.1840 0.5560];
    [0.4660 0.6740 0.1880];
    [0.3010 0.7450 0.9330]];
figure('Position', [548 458 692 338]), hold on
t = xt(rois_real, fs, 1);
for i = 1:size(rois_real,2)
    plot(t, zscore(rois_real(:,i))- i*5, 'k');
    plot(t, zscore(rois_modeled(:,i)) - i*5, 'Color', cols(i,:))
end
axis([15 315 -20 5])
axis off
line([51 60], [-19 -19], 'Color', [0 0 0], 'LineWidth', 2)
line([51 51], [-19 -17], 'Color', [0 0 0], 'LineWidth', 2)
% saveas(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'model_timeseries.svg']))
