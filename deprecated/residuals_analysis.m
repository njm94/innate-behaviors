% RESIDUALS ANALYSIS ON VIDEO RECONSTRUCTION

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

fs = 90 ;
[b, a] = butter(2, 0.01/(fs/2), 'high');

thy1_idx = 1:7;
ai94_idx = 8:13;
hyl3_idx = 14:19;
ibl2_idx = 22:25;

save_average_across_days = true;
num_comps = 200;
load('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat');
%%
 
for j = 1:length(data_list{1})+1
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
                                       
                savefig(gcf, [fPath,  char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_residuals.fig'])
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

        behavior_frames = cell(1, 17);
        
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
    events = read_boris(boris_file);

    % load experiment specific data into cell array
    load([data_dir filesep 'tform.mat'])

    disp('Loading brain data...')
    load(brain_file);
    U = permute(reshape(U, 128, 128, []), [2 1 3]);
    U = evaluate_tform(U, tformEstimate); % apply registration
    U = reshape(U, 128*128, []);
    U = U(:,1:num_comps);

    % light-OFF signal may be captured in the brain data, resulting in a
    % massive filter artifact - remove a few frames just in case
    trial_length = size(V, 2) - 5; 
    V = V(1:num_comps, 1:trial_length);

    % Brain video may have been resampled to wrong number of frames due to
    % inaccuracy in pre-processing step. Fix this by checking the length of
    % brain data and comparing to manually labeled VideoEnd event from
    % BORIS file. If there is a discrepancy, resample the brain data to the
    % length of the BORIS file to fix the issue.
    trial_length = find(events.("Video End"));
    if size(V,2) ~= trial_length
        V = resamplee(V', trial_length, size(V,2))';
    end

    Vbrain = recastV(Umaster(:,1:num_comps), U, s(1:num_comps, 1:num_comps), V(:, 1:trial_length));
    Vmaster = filtfilt(b, a, Vbrain')';

    Vfull = recastV(Umaster(1:num_comps, U, eye(num_comps, num_comps), Vfull(:, 1:trial_length)));

    % transform one last time to get into coordinates of atlas
%     U = permute(reshape(U, 128, 128, []), [2 1 3]);
%     U = evaluate_tform(U, atlas_tform.tform); % apply registration
%     U = reshape(U, 128*128, []);
    
    disp('Computing DF/F0')
    dataR = permute(reshape(Umaster*Vmaster, 128, 128, []), [2 1 3]);
    dFF = zscore((dataR-min(dataR(:))./mean(dataR - min(dataR(:)), 3))-1, [], 3);
    clear dataR Vbrain U s V

    % %% load grooming events
    % disp('Loading grooming events')
    % [behaviors, annotation_labels] = parse_snippets(snippets_dir);
    % 
    % left= idx2arr(behaviors{annotation_labels == "left"}, trial_length);
    % right= idx2arr(behaviors{annotation_labels == "right"}, trial_length);
    % elliptical= idx2arr(behaviors{annotation_labels == "elliptical"}, trial_length);
    % largeleft= idx2arr(behaviors{annotation_labels == "largeleft"}, trial_length);
    % largeright= idx2arr(behaviors{annotation_labels == "largeright"}, trial_length);
    % largebilateral= idx2arr(behaviors{annotation_labels == "largebilateral"}, trial_length);
    % lick = idx2arr(behaviors{annotation_labels == "lick"}, trial_length);


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

    end


    % consolidate lick events into a single variable in the table
    lick_idx = contains(events.Properties.VariableNames, 'Lick');

    Lick = events(:,lick_idx);
    events = removevars(events, lick_idx);
    Lick = any(table2array(Lick),2); 
    events = addvars(events, logical(Lick));



    bvars = ["Drop Hits Left", "Drop Hits Right", "Drop Hits Center", "Audio", ...
        "Lick", ...
        "Right", "Left", "Right Asymmetric", "Left Asymmetric" ...
        "Elliptical", "Elliptical Asymmetric", "Large Bilateral", ...
        "LeftMove", "RightMove"];

    behavior_frames = cell(1, length(bvars));

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
        behavior_frames{ii} = dFF(:,:,b_idx(1:end));
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

    save([data_dir, filesep, 'outputs', filesep, char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_behavior_frames.mat'], 'bvars', 'behavior_frames', '-v7.3')

    disp('Clearing variables')
    clear trials timestamps dFF largeright largeleft largebilateral lick left right LeftMove RightMove
%     if ~save_average_across_days
%         clear behavior_frames
%     end
    % count = count + 1;
%     if j == 3, break; end
%     master_SVD_file = [mouse_root_dir filesep 'masterSVD.mat'];
%     if ~isfile(master_SVD_file)



% %% load behavior video SVD

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

