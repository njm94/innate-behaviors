%% Spontaneous vs evoked average responses

%% average events

clc, clear, close all


addpath('/home/user/Documents/grooming/utils')


cluster_data = fix_path('Y:\nick\behavior\grooming\20241112141427_behavior_clustering.mat');
cluster_data = fix_path('Y:\nick\behavior\grooming\20241114092737_behavior_clustering.mat');
fileID = fopen('expt1_datalist.txt','r');

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = '';

fs = 90 ;
[b, a] = butter(2, 0.01/(fs/2), 'high');

thy1_idx = 1:7;
ai94_idx = 8:13;
hyl3_idx = 14:19;
ibl2_idx = 20:25;

save_average_across_days = true;
include_boris = true;
prune_outliers = false;

load(fix_path('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat'));

bvars = ["DropLeft", "DropRight", "DropCenter", "Audio", ...
        "Lick", ...
        "Right", "Left", "Right Asymmetric", "Left Asymmetric" ...
        "Elliptical", "Elliptical Left", "Elliptical Right" ...
        "LeftMove", "RightMove"];


%%
spon_index = [1 2 8 9 14 15 20 21];
evoked_index = 1:length(data_list{1});
evoked_index(spon_index) = [];

for j = 1:length(evoked_index)+1
     try

        data_dir = fix_path(data_list{1}{evoked_index(j)});
        disp(['Starting ' data_dir])

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
                disp('Saving average maps of individual behaviors...')
                save([fPath, char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), 'evoked_avg_behavior_frames.mat'], 'bvars', 'avg_behavior_frames', 'behavior_counter', '-v7.3')
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
                                       
                savefig(gcf, [fPath,  char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_evokeddFF.fig'])
            end
            %%
            % all mice completed - break the loop
            if j == length(evoked_index)+1, break; end
        end        
        count = 1;
        disp('Loading master basis set')
        load([mouse_root_dir filesep 'Umaster.mat'])
        load([mouse_root_dir filesep 'mask.mat'])
        atlas_tform = load([mouse_root_dir filesep 'atlas_tform.mat']);
        clear Vmaster left right elliptical large_left large_right largebilateral lick LeftMove RightMove Audio drop_left drop_right
        current_mouse = mouse_id;

        fPath = [mouse_root_dir filesep 'outputs' filesep];

        behavior_frames = cell(1, length(bvars));
        avg_behavior_frames = cell(1, length(bvars));
        behavior_counter = zeros(1, length(bvars));
        
    end

    try
        brain_file = [data_dir, filesep, get_file_with_str(data_dir, 'cam0_svd')];
        dlc_speed_file = [data_dir, filesep, get_file_with_str(data_dir, 'speed.csv')];
        timestamp_file = [data_dir, filesep, get_file_with_str(data_dir, 'trim.txt')];
        boris_file = [data_dir, filesep, get_file_with_str(data_dir, 'events.tsv')];
    catch
        disp('Missing one or more of the required files. Skipping...')
        continue
    end

    %% load grooming events from BORIS file
    [events, b_idx, ~, trial_length, cluster_labels] = get_labels_from_clustering_results(cluster_data, boris_file, include_boris, prune_outliers);

    % load experiment specific data into cell array
    load([data_dir filesep 'tform.mat'])


    disp('Loading brain data...')
    load(brain_file);
    U = permute(reshape(U, 128, 128, []), [2 1 3]);
    U = evaluate_tform(U, tformEstimate); % apply registration
    U = reshape(U, 128*128, []);

    % light-OFF signal may be captured in the brain data, resulting in a
    % massive filter artifact - remove a few frames just in case
    V = V(:,1:end-5);

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
    dataR = permute(reshape(Umaster*Vmaster, 128, 128, []), [2 1 3]);
    dFF = zscore((dataR-min(dataR(:))./mean(dataR - min(dataR(:)), 3))-1, [], 3);
    clear dataR Vbrain U s V

    %% load DLC tracks
    disp('Loading DLC tracks')
    dlc_speed = readmatrix(dlc_speed_file);
    fll_speed = dlc_speed(1:trial_length,1);
    flr_speed = dlc_speed(1:trial_length,2);
    
    % binarize movement speed
    LeftMove = fll_speed > mean(fll_speed) + std(fll_speed);
    RightMove = flr_speed >  mean(flr_speed) + std(flr_speed);

    % LeftMove = aggregate(fll_speed, )
    
    % k = 35; % smoothing kernel - this matches smoothing kernel in grooming detection but may need to change
    % LeftMove= movmax(fll_speed, k);
    % RightMove= movmax(flr_speed, k);

    clear dlc_speed fll_speed flr_speed

    %% load stimulus info from timestamp file
    % 10kHz tone occurs for 1 second immediately at start of trial
    % Tone is followed by 1 second of silence
    % Then valve is opened for 1.5s. Water drop occurs somewhere in there
    disp('Getting stimulus info from timestamp file');
    timestamps = readmatrix(timestamp_file);
    trials = unique(timestamps(:, 4));
    Audio= false(size(LeftMove));
    for i = 1:length(trials)
        if trials(i) == 0, continue; end
        trial_start  = find(timestamps(:,4)==trials(i), 1); 
        Audio(trial_start:trial_start+fs) = true;
    end

    % Iterate through behavior index vector. 
    % Extract behavior frames and put them in the behavior frame cell
    % Compute average of the frames within that behavior
    % Place that into the average behavior frame cell
    for ii = 1:size(b_idx,1)
        % figure out which cell to add the results in
        cell_idx = find(strcmpi(bvars, cluster_labels{ii}));

        behavior_frames{cell_idx} = cat(3, behavior_frames{cell_idx}, dFF(:,:,b_idx(ii,1):b_idx(ii,2)));
        behavior_counter(cell_idx) = behavior_counter(cell_idx)+1;
        avg_behavior_frames{cell_idx} = cat(3, avg_behavior_frames{cell_idx}, mean(dFF(:,:,b_idx(ii,1):b_idx(ii,2)),3));        
    end

    % place movement and audio into the cell array
    cell_idx = find(strcmpi(bvars, 'LeftMove'));
    behavior_frames{cell_idx} = cat(3, behavior_frames{cell_idx}, dFF(:,:,LeftMove));
    avg_behavior_frames{cell_idx} = cat(3, avg_behavior_frames{cell_idx}, mean(dFF(:,:,LeftMove),3));
    behavior_counter(cell_idx) = behavior_counter(cell_idx)+1;

    cell_idx = find(strcmpi(bvars, 'RightMove'));
    behavior_frames{cell_idx} = cat(3, behavior_frames{cell_idx}, dFF(:,:,RightMove));
    avg_behavior_frames{cell_idx} = cat(3, avg_behavior_frames{cell_idx}, mean(dFF(:,:,RightMove),3));
    behavior_counter(cell_idx) = behavior_counter(cell_idx)+1;

    cell_idx = find(strcmpi(bvars, 'Audio'));
    behavior_frames{cell_idx} = cat(3, behavior_frames{cell_idx}, dFF(:,:,Audio));
    avg_behavior_frames{cell_idx} = cat(3, avg_behavior_frames{cell_idx}, mean(dFF(:,:,Audio),3));
    behavior_counter(cell_idx) = behavior_counter(cell_idx)+1;

    % Place licks and drops into the cell array
    if any(contains(events.Properties.VariableNames, 'Lick'))
        lick_idx = arr2idx(events.Lick);
        cell_idx = find(strcmpi(bvars, 'Lick'));
        for ii = 1:size(lick_idx,1)
            behavior_frames{cell_idx} = cat(3, behavior_frames{cell_idx}, dFF(:,:,lick_idx(ii,1):lick_idx(ii,2)));
            avg_behavior_frames{cell_idx} = cat(3, avg_behavior_frames{cell_idx}, mean(dFF(:,:,lick_idx(ii,1):lick_idx(ii,2)),3));        
            behavior_counter(cell_idx) = behavior_counter(cell_idx)+1;
        end
    end

    if any(contains(events.Properties.VariableNames, 'Drop'))
        drop_idx = find(contains(events.Properties.VariableNames, 'Drop'));
        for ii = 1:length(drop_idx)
            cell_idx = find(strcmpi(bvars, events.Properties.VariableNames{drop_idx(ii)}));
            drop_event_idx = find(table2array(events(:,drop_idx(ii))));
            for jj = 1:size(drop_event_idx,1)
                % Drops are labelled at the time at which the full weight
                % of the drop lands on the face. This can mean that the
                % drop initially touches the face, but is still hanging on
                % the spout prior to the event. For this reason, take 1.5s
                % before the drop event, to 0.5s after.
                behavior_frames{cell_idx} = cat(3, behavior_frames{cell_idx}, dFF(:, :, drop_event_idx(jj)-round(1.5*fs):drop_event_idx(jj)+round(0.5*fs)));
                avg_behavior_frames{cell_idx} = cat(3, avg_behavior_frames{cell_idx}, mean(dFF(:, :, drop_event_idx(jj)-round(1.5*fs):drop_event_idx(jj)+round(0.5*fs)),3));
                behavior_counter(cell_idx) = behavior_counter(cell_idx)+1;
            end

        end
    end

    disp('Clearing variables')
    clear trials timestamps dFF largeright largeleft largebilateral lick left right LeftMove RightMove cell_idx

end


%% generate maps of spontaneous vs evoked
clc, clear, close all


addpath('/home/user/Documents/grooming/utils')

fileID = fopen('expt1_datalist.txt','r');

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = '';


thy1_idx = 1:7;
ai94_idx = 8:13;
hyl3_idx = 14:19;
ibl2_idx = 20:25;


load(fix_path('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat'));

bvars = ["DropLeft", "DropRight", "DropCenter", "Audio", ...
        "Lick", ...
        "Right", "Left", "Right Asymmetric", "Left Asymmetric" ...
        "Elliptical", "Elliptical Left", "Elliptical Right" ...
        "LeftMove", "RightMove"];


%%
spon_index = [1 2 8 9 14 15 20 21];
evoked_index = 1:length(data_list{1});
evoked_index(spon_index) = [];




%% spontaneous


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
    dff_fig = getAllFiles(dff_path, '_spontaneousdFF.fig');
    bFrames = getAllFiles(dff_path, 'spontaneous_avg');

    % If there are multiple versions, use the most recent one
    if size(dff_fig,1) > 1
        dff_fig = sort(dff_fig);
        dff_fig = dff_fig{end};
    end

    % If there are multiple versions, use the most recent one
    if size(bFrames,1) > 1
        bFrames = sort(bFrames);
        bFrames = bFrames{end};
    end
    load([dff_path, filesep, bFrames], 'behavior_counter');
    if j==1
        bcounter = behavior_counter;
    else
        bcounter = bcounter + behavior_counter;
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
bcounter = bcounter(5:end-2);
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

vars = ["lick", "right", "left", "largeright", "largeleft", "elliptical", "ellip_right", "ellip_left"];
% figure, axis off, hold on
for j = 1:length(mice)
    load(fix_path([data_root, filesep, mice{j}, filesep, 'atlas_tform.mat']));
    for i = 1:length(vars)
        
        behavior_map = eval(vars(i));
        if j > size(behavior_map,3)
            continue
        end
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
%                 title([char(vars(i)), ' n=', num2str(bcounter(i))]);
            else
%                 disp([char(vars(i)), ' is not greater than 0'])
            end
            set(gca, 'YDir', 'reverse');
    end
    
end
% % 
% for i = 1:length(vars)
%     disp([char(vars(i)), ' n=', num2str(bcounter(i))]);
%     exportgraphics(figure(i), fix_path(['Y:\nick\behavior\grooming\figures\',char(vars(i)), '_contours_spontaneous.png']), 'Resolution', 300)
% end

%% evoked

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
    dff_fig = getAllFiles(dff_path, '_evokeddFF.fig');
    bFrames = getAllFiles(dff_path, 'evoked_avg');

    % If there are multiple versions, use the most recent one
    if size(dff_fig,1) > 1
        dff_fig = sort(dff_fig);
        dff_fig = dff_fig{end};
    end

    % If there are multiple versions, use the most recent one
    if size(bFrames,1) > 1
        bFrames = sort(bFrames);
        bFrames = bFrames{end};
    end
    load([dff_path, filesep, bFrames], 'behavior_counter');
    if j==1
        bcounter = behavior_counter;
    else
        bcounter = bcounter + behavior_counter;
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
bcounter = bcounter(5:end-2);
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

vars = ["lick", "right", "left", "largeright", "largeleft", "elliptical", "ellip_right", "ellip_left"];
% figure, axis off, hold on
for j = 1:length(mice)
    load(fix_path([data_root, filesep, mice{j}, filesep, 'atlas_tform.mat']));
    for i = 1:length(vars)
        
        behavior_map = eval(vars(i));
        if j > size(behavior_map,3)
            continue
        end
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
%                 title([char(vars(i)), ' n=', num2str(bcounter(i))]);
            else
%                 disp([char(vars(i)), ' is not greater than 0'])
            end
            set(gca, 'YDir', 'reverse');
    end
    
end
%% 
% for i = 1:length(vars)
%     disp([char(vars(i)), ' n=', num2str(bcounter(i))]);
%     exportgraphics(figure(i), fix_path(['Y:\nick\behavior\grooming\figures\',char(vars(i)), '_contours_evoked.png']), 'Resolution', 300)
% end
