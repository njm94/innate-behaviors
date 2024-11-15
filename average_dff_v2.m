%% Make average responses

%% average events

clc, clear


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
 
for j = 1:length(data_list{1})+1
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
        atlas_tform = load([mouse_root_dir filesep 'atlas_tform.mat']);
        clear Vmaster left right elliptical large_left large_right largebilateral lick LeftMove RightMove Audio drop_left drop_right
        current_mouse = mouse_id;

        fPath = [mouse_root_dir filesep 'outputs' filesep];

        behavior_frames = cell(1, length(bvars));
        
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
    [events, ~, ~, trial_length] = get_labels_from_clustering_results(cluster_data, boris_file, include_boris, prune_outliers);

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
    fll_speed = fll_speed > mean(fll_speed) + std(fll_speed);
    flr_speed = flr_speed >  mean(flr_speed) + std(flr_speed);
    
    k = 35; % smoothing kernel - this matches smoothing kernel in grooming detection but may need to change
    LeftMove= movmax(fll_speed, k);
    RightMove= movmax(flr_speed, k);

    clear dlc_speed fll_speed flr_speed

    %% load stimulus info from timestamp file
    % 10kHz tone occurs for 1 second immediately at start of trial
    % Tone is followed by 1 second of silence
    % Then valve is opened for 1.5s. Water drop occurs somewhere in there
    disp('Getting stimulus info from timestamp file');
    timestamps = readmatrix(timestamp_file);
    trials = unique(timestamps(:, 4));
    Audio= zeros(size(LeftMove));
    % drop_left = zeros(size(fll_move));
    % drop_right = zeros(size(fll_move));
    for i = 1:length(trials)
        if trials(i) == 0, continue; end
        trial_start  = find(timestamps(:,4)==trials(i), 1); 
        Audio(trial_start:trial_start+fs) = 1;
        % if count > 2 && count < 6
        %     % first 2 trials for each mouse are spontaneous trials - no
        %     % water drop. then on last trial for HYL3, drop goes on left
        %     % side
        %     drop_right(trial_start+ceil(2*fs): trial_start+ceil(3.5*fs))=1;
        % elseif count == 6
        %     drop_left(trial_start+ceil(2*fs): trial_start+ceil(3.5*fs))=1;
        % end
    end

    % iterate through behaviors and store frames
    

    % Update drop variable to capture window surrounding each event
    for ii = 1:3
        if any(strcmp(events.Properties.VariableNames, bvars(ii)))
            drop_event_idx = find(table2array(events(:,ii)));
            for jj = 1:length(drop_event_idx)
                events.(bvars(ii))(drop_event_idx(jj)-30:drop_event_idx(jj)+90) = 1;
            end
        end
    end

    figure
    for ii = 1:length(bvars)
        if any(strcmp(events.Properties.VariableNames, bvars(ii)))
            b_idx = logical(aggregate(table2array(events(:,strcmp(events.Properties.VariableNames, bvars(ii)))), 3));
        elseif exist(bvars(ii), 'var')
            if any(eval(bvars(ii)))
                b_idx = logical(aggregate(eval(bvars(ii)),3));
            else
                continue
            end
        else
            continue
        end
        behavior_frames{ii} = cat(3, behavior_frames{ii}, dFF(:,:,b_idx));
        subplot(5, 4, ii)

        % transpose the mean image to get the brains in proper
        % orientation and % transform to get into coordinates of atlas
        mean_image = mask.*mean(dFF(:,:,b_idx(1:end)),3)';
        mean_image = imwarp(mean_image, atlas_tform.tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
        imagesc(mean_image);
        title(bvars(ii)), colorbar, xticks([]),  yticks([])
        hold on;
        for p = 1:length(dorsalMaps.edgeOutline)
            plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'w');
        end
        set(gca, 'YDir', 'reverse');
    end

    savefig(gcf, [data_dir, filesep, 'outputs', filesep,  char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_dFF.fig'])
    close(gcf)

    % disp('Saving...')
    % save([data_dir, filesep, 'outputs', filesep, char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_behavior_frames.mat'], 'bvars', 'behavior_frames', '-v7.3')

    disp('Clearing variables')
    clear trials timestamps dFF largeright largeleft largebilateral lick left right LeftMove RightMove

end


