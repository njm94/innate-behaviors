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
else
    addpath('/home/user/Documents/grooming/ridgeModel');
    addpath('/home/user/Documents/grooming/ridgeModel\widefield')
    addpath('/home/user/Documents/grooming/ridgeModel\smallStuff') 
    addpath('/home/user/Documents/grooming/utils')
end




fileID = fopen('expt1_datalist.txt','r');

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = '';

fs = 90 ; % sampling rate
[b, a] = butter(2, 0.01/(fs/2), 'high');
aggregation_sz = 3; % window of time to aggregate behaviors
clen = 15; % duration after which continuous behavior is considered long
blen = 10; % duration of baseline period before and after continuous grooming event

thy1_idx = 1:7;
ai94_idx = 8:13;
camk_idx = 14:25;

save_average_across_days = true;

load('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat');
%%
 
for j = thy1_idx%camk_idx %23:length(data_list{1})+1
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

    if ~any(episode_durations >= 15)
%         disp(['No episodes longer than ', num2str(clen),'s. Skipping...'])
        continue
    else
        idx = idx(episode_durations >= 15, :);
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

%         tmp_idx(i,:) = [idx(i,1)-round(blen*fs), idx(i,2)+round(blen*fs)];
        dFF_crop = dFF(:,:,idx(1):idx(2));

%         dFF_crop = dFF_crop.*mask;
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
        corrmat{j}(:,:,i) = corrcoef(ts);

%         % perform a moving window correlation over the data
%         clear corrmat
%         stepsize = round(0.25 * fs);
%         corrk = 1 * fs;
%         count = 1;
%         for ii = 1:stepsize:size(ts,1)-corrk
%             corrmat(:,:,count) = corrcoef(ts(ii:ii+corrk,:));
%             count = count +1;
%         end

        % find a period of time that is equal in length with no activity
        groom_dur = size(ts,1);
        
        a = strfind(tmp', zeros(1,groom_dur));
        if isempty(a)
            continue
        else
            dFF_crop = dFF(:,:,a(1):a(1)+groom_dur);
            dFF_crop = imwarp(dFF_crop, atlas_tform.tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));

            ts_rest = getTimeseries(dFF_crop, seeds, 2);
            rest_corrmat{j}(:,:,i) = corrcoef(ts_rest);

            tmp(a(1):a(1)+groom_dur)=1;
            count = count + 1;
        end

    end


end


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