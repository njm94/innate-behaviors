%% Make average responses

%% average events


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



clc, clear
fileID = fopen('expt1_datalist.txt','r');

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = '';

fs = 90 ;
[b, a] = butter(2, 0.01/(fs/2), 'high');

thy1_idx = 1:7;
ai94_idx = 8:13;
camk_idx = 14:25;

save_average_across_days = false;
%%
 
for j = thy1_idx(end) %23:length(data_list{1})+1
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

            if save_average_across_days
                figure
                for iii = 1:size(behavior_frames,2)
                    % transpose the mean image to get the brains in proper
                    % orientation
                    mean_image = mean(behavior_frames{iii},3)';
                    subplot(3,4,iii)
                    imagesc(mean_image)
                    set(gca, 'clim', [-max(abs(mean_image(:))) max(abs(mean_image(:)))])
                    c=colorbar;
                    c.Label.String = '\DeltaF/F_0 (\sigma)';
                    title(bvars(iii))
                    xticks([])
                    yticks([])
                    colormap(fliplr(redblue([],[],'k')))
                end
                                       
                savefig(gcf, [fPath,  char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_dFF.fig'])
            end
            % all mice completed - break the loop
            if j == length(data_list{1})+1, break; end
        end        
        count = 1;
        disp('Loading master basis set')
        load([mouse_root_dir filesep 'Umaster.mat'])
        load([mouse_root_dir filesep 'mask.mat'])
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

    % load experiment specific data into cell array
    load([data_dir filesep 'tform.mat'])


    disp('Loading brain data...')
    load(brain_file);
    U = permute(reshape(U, 128, 128, []), [2 1 3]);
    U = evaluate_tform(U, tformEstimate); % apply registration
    U = reshape(U, 128*128, []);

    % Brain video may have been resampled to wrong number of frames due to
    % inaccuracy in pre-processing step. Fix this by checking the length of
    % brain data and comparing to manually labeled VideoEnd event from
    % BORIS file. If there is a discrepancy, resample the brain data to the
    % length of the BORIS file to fix the issue.
    trial_length = find(events.("Video End"));
    if size(V,2) > trial_length
        V = resamplee(V', trial_length, size(V,2))';
    end

    
    % light-OFF signal may be captured in the brain data, resulting in a
    % massive filter artifact - remove a few frames just in case
    trial_length = size(V, 2) - 5; 
    Vbrain = recastV(Umaster, U, s, V(:, 1:trial_length));
    Vmaster = filtfilt(b, a, Vbrain')';
    
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
    % bvars = ["elliptical", "largeleft", "largeright", "largebilateral", ...
    %     "left", "right", "lick", "fll_move", "flr_move", ...
    %     "audio_tone", "drop_left", "drop_right"];

    bvars = ["Drop Hits Left", "Drop Hits Right", "Drop Hits Center", "Audio", ...
        "Lick Left", "Lick Right", "Lick Unknown", "Lick No Paws", ...
        "Right", "Left", "Right Asymmetric", "Left Asymmetric" ...
        "Elliptical", "Elliptical Asymmetric", "Large Bilateral", ...
        "LeftMove", "RightMove"];

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
            b_idx = logical(table2array(events(:,strcmp(events.Properties.VariableNames, bvars(ii)))));
        elseif exist(bvars(ii), 'var')
            b_idx = logical(eval(bvars(ii)));
        else
            continue
        end
        behavior_frames{ii} = cat(3, behavior_frames{ii}, dFF(:,:,b_idx(1:end-5)));
        subplot(5, 4, ii)
        % transpose the mean image to get the brains in proper
        % orientation
        imagesc(mask.*mean(dFF(:,:,b_idx(1:end-5)),3)');
        title(bvars(ii)), colorbar, xticks([]),  yticks([])
    end

    % figure,  
    % for ii = 1:length(bvars)
    %     behavior_frames{ii} = cat(3, behavior_frames{ii}, dFF(:,:,logical(eval(bvars(ii)))));
    %     if any(logical(eval(bvars(ii))))
    %         subplot(4, 3, ii)
    %         % transpose the mean image to get the brains in proper
    %         % orientation
    %         imagesc(mean(dFF(:,:,logical(eval(bvars(ii)))), 3)');
    %         title(bvars(ii)), colorbar, xticks([]), yticks([])
    %     end
    % end
    savefig(gcf, [data_dir, filesep, 'outputs', filesep,  char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_dFF.fig'])
    close(gcf)

    disp('Clearing variabls')
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