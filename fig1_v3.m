% This code will read the clustered behaviors and generate figure 1
addpath('/home/user/Documents/grooming/utils')
load(fix_path('Y:\nick\behavior\grooming\20241114092737_behavior_clustering.mat'))

% map dendrogram labels to behaviors - Note this is subject to change based
% on the ordering of the output labels from the dendrogram in python
% script. Double check this.
label_map = ["Right", "Left", "Left Asymmetric", ...
    "Elliptical Asymmetric", "Right Asymmetric", "Elliptical"];

figure, hold on, 
% show UMAP embedding

labs = unique(dendrogram_labels);
for i = 1:length(labs)
    scatter(embedding(dendrogram_labels==labs(i),1), embedding(dendrogram_labels==labs(i),2),'.')
    axis off equal
end
ax = gca;

exportgraphics(ax, fix_path(['Y:\nick\behavior\grooming\figures\','UMAP', '.png']), 'Resolution', 300)
saveas(ax, fix_path(['Y:\nick\behavior\grooming\figures\','UMAP', '.svg']))

%%

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
    data_root = fix_path(fileparts(behavior_files{i,:}));
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
            if counter(j)>50, continue; end
            plot(nose_ref_x(j)-151-flr_x(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), nose_ref_y(j)+flr_y(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), 'Color', [1 0 1 0.15])
            plot(nose_ref_x(j)-151-fll_x(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), nose_ref_y(j)+fll_y(bIdx(clus_j(k),1):bIdx(clus_j(k),2)), 'Color', [0 1 1 0.15])
        end
    end
end


%%

for i = 1:length(label_map)
    figure(i);
    % axis([0, 300, 0, 300])
    ax = gca;
%     exportgraphics(ax, append('Y:\nick\behavior\grooming\figures', label_map(i),'.eps'), 'ContentType', 'vector')
    % exportgraphics(ax, append('Y:\nick\behavior\grooming\figures\', label_map(i),'.png'), 'Resolution', 300)
    % saveas(ax, fix_path(append('Y:\nick\behavior\grooming\figures\',label_map(i), '.pdf')))

end

%% get data for behavior event quantifications

fileID = fopen('expt1_datalist.txt','r');

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = ''; 
% addpath('C:\Users\user\Documents\Nick\grooming\deprecated')

% Notes on 2p data
% ETR2 and ETR3 had corrupted videos for spontaneous trials, so we will
% just use evoked videos
%
% ECL3, IDR3, and RR3 have 3 days of spontaneous behavior - use last 2 days
% for the timeline plot. The spontaneous data may not have Boris labelling.
% Use the snippets data instead. Also these mice have longer baselines in
% the evoked trials. (5mins in these mice vs 30s in the 1p mice). 
% They don't do any grooming anyway so should be comparable
%
% ECL3 20240729 corrupted synchronization between brain and behavior. Can't
% use the brain data with grooming together but grooming alone should be ok

mp_list = {'Y:\nick\behavior\grooming\2p\ETR2_thy1\20231113143925'; 
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231113155903'; 
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231115174148'; 
    '/media/user/teamshare/nick/behavior/grooming/2p/ECL3_thy1/20240722/';
    '/media/user/teamshare/nick/behavior/grooming/2p/ECL3_thy1/20240724/';
    '/media/user/teamshare/nick/behavior/grooming/2p/ECL3_thy1/20240726/';
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240729'; 
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240731';
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240802';
    '/media/user/teamshare/nick/behavior/grooming/2p/IDR3_tTA6s/20240722/';
    '/media/user/teamshare/nick/behavior/grooming/2p/IDR3_tTA6s/20240724/';
    '/media/user/teamshare/nick/behavior/grooming/2p/IDR3_tTA6s/20240726/';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240729';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240731';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240802';
    '/media/user/teamshare/nick/behavior/grooming/2p/RR3_tTA8s/20240722/';
    '/media/user/teamshare/nick/behavior/grooming/2p/RR3_tTA8s/20240724/';
    '/media/user/teamshare/nick/behavior/grooming/2p/RR3_tTA8s/20240726/';
    'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240729';
    'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240802'};

data_list{1} = [data_list{1}; mp_list];
for i = 1:size(data_list{1},1)
    data_list{1}{i} = fix_path(data_list{1}{i});
end


% make behavior event raster
spon_index = [1,2,8,9,14,15,20,21,29:31,35:37,41:43];
evoked_index = ones(size(data_list{1}));
evoked_index(spon_index) = 0;
evoked_index = find(evoked_index);

% % The following videos had issues upon acquisition where they are about 1
% % minute shorter than the rest. This will make it difficult to visualize
% % them in the raster since they are not properly aligned. Exclude from
% % visualization, but we will keep them in the quanitification
% exclude_for_raster = [26, 32, 42];

%%

aggregation_sz = 3;

N = length(data_list{1});

num_episodes = zeros(N, 1);


% use the Boris labels (old script) since we don't care about the type of
% grooming behavior here - just interested in episodes
for j = 1:N
    data_dir = data_list{1}{j};
    timestamp_file = [data_dir, filesep, getAllFiles(data_dir, '_trim.txt')];
    [mouse_root_dir, exp_date, ~] = fileparts(data_dir);
    [~, mouse_id, ~] = fileparts(mouse_root_dir);
    [expt_root_dir, ~, ~] = fileparts(mouse_root_dir);
    snippets_dir = [data_dir, filesep, 'snippets'];

    timestamps = readmatrix(timestamp_file);
    first_trial(j) = find(timestamps(:,3),1);

    if ~isempty(getAllFiles(data_dir, '.tsv'))
        boris_file = [data_dir, filesep, getAllFiles(data_dir, '.tsv')];
        [events, snippets, b_table] = read_boris(boris_file);

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
        bmat = any(table2array(stroke_events),2);
        % snippets = arr2idx(bmat);
    else 
        
        [snippets, labels] = parse_snippets(snippets_dir);

        % remove licks
        snippets(strcmpi(labels, 'lick')) = [];
        bmat = idx2arr(catcell(1, snippets));

    end

    % num_events(j) = cellfun(@(data) size(data, 1), snippets);
    if isempty(bmat)
        event_raster{j} = 0;
        num_episodes(j) = 0;
    else
        [event_raster{j}, idx] = aggregate(bmat, aggregation_sz);
        num_episodes(j) = size(idx,1);
    end
    
end


% %% align all trials
% 

% % that should all remain the same
% tlen = cellfun(@(x) size(x, 1), event_raster);
% [tlen, I] = median(tlen);
% 
% 
% clear raster_mat
% 
% % align everything to latest first trial
% tstart = max(first_trial);
% 
% raster_mat = zeros(length(data_list{1})-length(exclude_for_raster), tlen);
% count = 1;
% for i = 1:length(spon_index)
%     if any(exclude_for_raster==i)
%         continue
%     end
%     iStart = first_trial(spon_index(i));
%     iRaster = event_raster{spon_index(i)};
%     if length(iRaster) == 1
%         raster_mat(count, :) = 0;
%     else
%         raster_mat(count, tstart-(iStart-1):tstart+length(iRaster)-iStart) = iRaster;
%     end
%     count = count + 1;
% end
% figure, imagesc(1-raster_mat), colormap gray
% %%
% for i = 1:length(evoked_index)
%     if any(exclude_for_raster==i)
%         continue
%     end
%     iStart = first_trial(evoked_index(i));
%     iRaster = event_raster{evoked_index(i)};
%     if length(iRaster) == 1
%         raster_mat(count, :) = 0;
%     else
%         raster_mat(count, tstart-(iStart-1):tstart+length(iRaster)-iStart) = iRaster;
%     end
%     count = count + 1;
% end
% figure, imagesc(1-raster_mat), colormap gray

%% Align all trials

% the last cohort of 2p data mice have 5 minute baseline periods rather
% than 30s baseline periods like all the other mice. Use the 30s period, so
% the raster is comparable across all animals. Align them by the end, and
% they will line up better

tlen = cellfun(@(x) size(x, 1), event_raster);
tlen = median(tlen);


clear raster_mat

for i  = 1:length(event_raster)
    if length(event_raster{i}) == 1 % no grooming instances in the session
        raster_mat(i,:) = zeros(1, tlen);
    elseif length(event_raster{i}) < tlen
        raster_mat(i,:) = [zeros(tlen-length(event_raster{i}),1); event_raster{i}];
    else
        raster_mat(i,:) = event_raster{i}(end-tlen+1:end);
    end
end


%%
close all, clc
figure%('Position', [380 458 860 406]),
% subplot(2,1,1)
imagesc(1-[raster_mat(spon_index,:); raster_mat(evoked_index,:)])
colormap gray, 
axis off
y = [length(spon_index)+0.5, length(spon_index)+0.5, size(raster_mat,1)+0.5, size(raster_mat,1)+0.5];
patch([0 size(raster_mat,2) size(raster_mat,2) 0], y, [0.3010 0.7450 0.9330], 'FaceAlpha', 0.1, 'EdgeColor', 'none')
% hold on

ax = gca;
% set(gcf, 'PaperPositionMode', 'auto')
% set(ax,'renderer','painters');

% exportgraphics(ax, fix_path(['Y:\nick\behavior\grooming\figures\', 'raster', '.png']), 'Resolution', 300)
% saveas(ax, fix_path(['Y:\nick\behavior\grooming\figures\','raster1', '.svg']))
%%
p_evoked = smoothdata(mean(raster_mat(evoked_index,:)));
p_spon = smoothdata(mean(raster_mat(spon_index,:)));
t = xt(p_evoked, fs);
figure, hold on
plot(t, p_evoked, 'LineWidth', 2)
plot(t, p_spon, 'k', 'LineWidth', 2)

% drop = 30:60:t(end);
% for i = 1:length(drop)
%     vline(drop, 'k:')
% end
axis tight
xlabel('Time (s)')
ylabel('Probability')

ax = gca;
exportgraphics(ax, fix_path(['Y:\nick\behavior\grooming\figures\', 'psth', '.png']), 'Resolution', 300)
saveas(ax, fix_path(['Y:\nick\behavior\grooming\figures\','psth', '.svg']))
%%

figure
num_spon_events = num_episodes(spon_index);
num_evoked_events = num_episodes(evoked_index);
[h,p] = ttest2(num_spon_events, num_evoked_events);

size_diff = length(num_evoked_events) - length(num_spon_events);
num_spon_events = cat(1, num_spon_events, nan(size_diff, 1));
data2plot = [num_spon_events, num_evoked_events];

boxplot(data2plot, 'colors', 'kk', 'Labels',{'Spontaneous','Evoked'}, 'Symbol', '');
ylabel('# Grooming episodes')
hold on, swarmchart(ones(size(num_spon_events))+0.3, num_spon_events, 'ko', 'XJitterWidth', 0.25)
swarmchart(ones(size(num_evoked_events))+1.3, num_evoked_events, 'ko', 'XJitterWidth', 0.25)
ax = gca;
ax.FontSize = 12;
title(['p=', num2str(p)])

% ax = gca;
% exportgraphics(ax, fix_path(['Y:\nick\behavior\grooming\figures\', 'num_events_sponvsevoked', '.png']), 'Resolution', 300)
% saveas(ax, fix_path(['Y:\nick\behavior\grooming\figures\','num_events_sponvsevoked', '.svg']))

%%

grouped_label_map = ["Elliptical", "Asymmetric", "Unilateral", "Elliptical Asymmetric"];

grouped_labels = zeros(size(dendrogram_labels));
grouped_labels(dendrogram_labels==find(strcmp(label_map, 'Elliptical'))) = find(strcmp(grouped_label_map, 'Elliptical'));
grouped_labels(dendrogram_labels==find(strcmp(label_map, 'Right Asymmetric')) | dendrogram_labels==find(strcmp(label_map, 'Left Asymmetric'))) = find(strcmp(grouped_label_map, 'Asymmetric'));
grouped_labels(dendrogram_labels==find(strcmp(label_map, 'Right')) | dendrogram_labels==find(strcmp(label_map, 'Left'))) = find(strcmp(grouped_label_map, 'Unilateral'));
grouped_labels(dendrogram_labels==find(strcmp(label_map, 'Elliptical Asymmetric'))) = find(strcmp(grouped_label_map, 'Elliptical Asymmetric'));



all_spon_index = zeros(size(dendrogram_labels));
bFileCell = cellstr(bFiles);
for i = 1:length(data_list{1})
    if any(spon_index == i)
        all_spon_index(contains(bFileCell, data_list{1}{i})) = 1;
    end   
end
all_evoked_index = 1-all_spon_index;


spon_stats = tabulate(grouped_labels(logical(all_spon_index)));
evoked_stats = tabulate(grouped_labels(logical(all_evoked_index)));

figure('Position', [675 453 293 413]), 
b=bar([spon_stats(:,3) evoked_stats(:,3)]);
b(1).FaceColor = 'k';
b(2).FaceColor =  [0.3010 0.7450 0.9330];

hold on
% errorbar((1:4)-0.15, mean(all_spon_relfreq,2), std(all_spon_relfreq,[],2)./sqrt(size(all_spon_relfreq,2)), 'k.')
% errorbar((1:4)+0.15, mean(all_evoked_relfreq,2), std(all_evoked_relfreq,[],2)./sqrt(size(all_evoked_relfreq,2)), 'k.')
xticklabels(grouped_label_map)
ylabel('Relative Frequency')
legend(["Spontaneous", "Evoked"], 'Location', 'NorthEast', 'Box', 'off')
ylim([0 100])
ax = gca;
ax.FontSize = 12;

ax = gca;
exportgraphics(ax, fix_path(['Y:\nick\behavior\grooming\figures\', 'rel_freq', '.png']), 'Resolution', 300)
saveas(ax, fix_path(['Y:\nick\behavior\grooming\figures\','rel_freq', '.svg']))