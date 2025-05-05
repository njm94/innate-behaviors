
clear, close all
cd('/media/user/teamshare/nick/behavior/grooming/code/')
addpath('/media/user/teamshare/nick/behavior/grooming/code/ridgeModel')
addpath('/media/user/teamshare/nick/behavior/grooming/code/ridgeModel/widefield')
addpath('/media/user/teamshare/nick/behavior/grooming/code/ridgeModel/smallStuff')

% addpath('C:\Users\user\Documents\Nick\ridgeModel');
% addpath('C:\Users\user\Documents\Nick\ridgeModel\widefield')
% addpath('C:\Users\user\Documents\Nick\ridgeModel\smallStuff') 

%%
clc
fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = '';

fs = 90;
[b, a] = butter(2, 0.01/(fs/2), 'high');

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
            % if all the trials have completed for the previous mouse, run
            % ridge regression
                       
            % select minimum length trial and then concatenate to array
            disp('Truncating trials to same length')
            min_trial_length = min(cellfun('size', Vmaster, 2));
            num_trials = length(Vmaster);
            for i = 1:num_trials
                Vmaster{i} = Vmaster{i}(:,1:min_trial_length);
                ipsi{i} = ipsi{i}(1:min_trial_length);
                contra{i} = contra{i}(1:min_trial_length);
                bilat{i} = bilat{i}(1:min_trial_length);
                fll_move{i} = fll_move{i}(1:min_trial_length);
                flr_move{i} = flr_move{i}(1:min_trial_length);
                audio_tone{i} = audio_tone{i}(1:min_trial_length);
                water_drop{i} = water_drop{i}(1:min_trial_length);
            end
            Vmaster = catcell(2, Vmaster);
            ipsi = catcell(2, ipsi);
            contra = catcell(2, contra);
            bilat = catcell(2, bilat);
            fll_move = catcell(1, fll_move)';
            flr_move = catcell(1, flr_move)';
            audio_tone = catcell(1, audio_tone)';
            water_drop = catcell(1, water_drop)';
            new_trial = repmat([1 zeros(1, min_trial_length-1)], 1, num_trials);

            % ridge here
            disp('Building design matrix')
            opts.frameRate = fs;
            opts.sPostTime=round(fs*2);
            opts.mPreTime = ceil(0.5 * fs);  % precede motor events to capture preparatory activity in frames (used for eventType 3)
            opts.mPostTime = ceil(2 * fs);   % follow motor events for mPostStim in frames (used for eventType 3)
            opts.framesPerTrial = min_trial_length; % nr. of frames per trial
            opts.folds = 10; %nr of folds for cross-validation
            
            % regressor_mat = [new_trial' audio_tone' water_drop' ipsi' contra' bilat' fll_move' flr_move']; % tone drop ipsi contra bilatInclusive forelimbInc
            regressor_mat = [audio_tone' water_drop' ipsi' contra' bilat' fll_move' flr_move']; % tone drop ipsi contra bilatInclusive forelimbInc
            % Full-Trial events:    new_trial
            % Post-Stimulus events: audio_tone, water_drop
            % Peri-Stimulus events: ipsi, contra, bilat, fll_move, flr_move
            % [dMat, regIdx] = makeDesignMatrix(regressor_mat, [1, 2, 2, 3, 3, 3, 3, 3], opts);
            [dMat, regIdx] = makeDesignMatrix(regressor_mat, [2, 2, 3, 3, 3, 3, 3], opts);
            regLabels = {'audio', 'drop', 'ipsi', 'contra', 'bilat', 'left_move', 'right_move'}; %some movement variables
            fullR = [dMat];

            disp('Running ridge regression with 10-fold cross-validation')
            [Vfull, fullBeta, ~, fullIdx, fullRidge, fullLabels] = crossValModel(fullR, Vmaster, regLabels, regIdx, regLabels, opts.folds);
            save([fPath 'cvFull.mat'], 'Vfull', 'fullBeta', 'fullR', 'fullIdx', 'fullRidge', 'fullLabels', '-v7.3'); %save some results
            
            fullMat = modelCorr(Vmaster,Vfull,Umaster) .^2;

            disp('Running reduced models')
            reducedMat = [];
            for i = 1:length(regLabels)
                reduced = fullR;
                cIdx = regIdx == i;
                if ~any(cIdx)
                    disp(['No ', regLabels{i}, '  events found. Skipping...'])
                    continue
                end
                reduced(:, cIdx) = reduced(randperm(size(reduced, 1)), cIdx);
            
                [Vreduced{i}, reducedBeta{i}, reducedR, reducedIdx, reducedRidge, reducedLabels] = crossValModel(reduced, Vmaster, regLabels, regIdx, regLabels, opts.folds);
                reducedMat(:, i) = modelCorr(Vmaster, Vreduced{i}, Umaster) .^2; %compute explained variance
            end

            save([fPath 'cvReduced.mat'], 'Vreduced', 'reducedBeta', 'reducedR', 'reducedIdx', 'reducedRidge', 'reducedLabels', '-v7.3'); %save some results

            figure
            subplot(3, 3, 1)
            imagesc(reshape(fullMat, [128 128]))
            xticks([])
            c=colorbar;
            title('FullMat')
            c.Label.String = 'cvR^2';
            yticks([])
            for i = 1:length(regLabels)
                subplot(3, 3, i+1)
                imagesc(reshape(fullMat - reducedMat(:,i), [128 128]))
                c=colorbar;
                title(regLabels{i})
                c.Label.String = '\DeltaR^2';
            
                xticks([])
                yticks([])
            end

        
            savefig(gcf, [fPath 'summary.fig'])


            % all mice completed - break the loop
            if j == length(data_list{1})+1, break; end
        end        
        count = 1;
        disp('Loading master basis set')
        load([mouse_root_dir filesep 'Umaster.mat'])
        clear Vmaster ipsi contra bilat fll_move flr_move audio_tone water_drop
        current_mouse = mouse_id;

        %     %% draw mask - might be needed for memory management
        %     frame = loadtiff(frame_file);
        %     mask = draw_roi(frame, 4);
    else
        fPath = [mouse_root_dir filesep 'outputs' filesep];
    end

    try
        brain_file = [data_dir, filesep, get_file_with_str(data_dir, 'cam0_svd')];
%         beh_file = [data_dir, filesep, get_file_with_str(data_dir, [exp_date, '_svd'])];
%         ME_file = [data_dir, filesep, get_file_with_str(data_dir, 'MEsvd')];
%         frame_file = [data_dir, filesep, get_file_with_str(data_dir, 'singleFrame')];
        grooming_file = [data_dir, filesep, 'grooming_events_roi_filtered.mat'];
        dlc_speed_file = [data_dir, filesep, get_file_with_str(data_dir, 'speed.csv')];
        timestamp_file = [data_dir, filesep, get_file_with_str(data_dir, 'trim.txt')];
    catch
        disp('Missing one or more of the required files. Skipping...')
        continue
    end
    
    % load experiment specific data into cell array
    load([data_dir filesep 'tform.mat'])

    disp('Loading brain data...')
    load(brain_file);
    U = permute(reshape(U, 128, 128, []), [2 1 3]);
    U = evaluate_tform(U, tformEstimate); % apply registration
    U = reshape(U, 128*128, []);

    Vbrain = recastV(Umaster, U, s, V);
    % Vmaster{count} = Vbrain;
    Vmaster{count} = filtfilt(b, a, Vbrain')';
    clear U s V Vbrain

    %% load grooming events
    disp('Loading grooming events')
    load(grooming_file)

    %     % don't treat bilateral as a separate event
%     ipsi = FL_R;
%     contra = FL_L;
    
    % separate ipsilateral, contralateral, and bilateral events
    ipsi{count} = zeros(size(FL_L));
    contra{count} = zeros(size(FL_L));
    bilat{count} = zeros(size(FL_L));
    for i = 1:size(event_type, 1)
        if contains(event_type(i,:), 'ipsi')
            ipsi{count}(all_events(i,1):all_events(i,2)) = 1;
        elseif contains(event_type(i,:), 'cont')
            contra{count}(all_events(i,1):all_events(i,2)) = 1;
        elseif contains(event_type(i,:), 'bi')
            bilat{count}(all_events(i,1):all_events(i,2)) = 1;
        else
            disp('wrong string in event type')
        end
    end

    %% load DLC tracks
    disp('Loading DLC tracks')
    dlc_speed = readmatrix(dlc_speed_file);
    fll_speed = dlc_speed(:,1);
    flr_speed = dlc_speed(:,2);
    
    % binarize movement speed
    fll_speed = fll_speed > mean(fll_speed) + std(fll_speed);
    flr_speed = flr_speed >  mean(flr_speed) + std(flr_speed);
    
    k = 35; % smoothing kernel - this matches smoothing kernel in grooming detection but may need to change
    fll_move{count} = movmax(fll_speed, k);
    flr_move{count} = movmax(flr_speed, k);

%     % eliminate overlap of movement and grooming
%     fll_move_ex = fll_move & ~(ipsi | contra | bilat)';
%     flr_move_ex = flr_move & ~(ipsi | contra | bilat)';

    clear dlc_speed fll_speed flr_speed

    %% load stimulus info from timestamp file
    % 10kHz tone occurs for 1 second immediately at start of trial
    % Tone is followed by 1 second of silence
    % Then valve is opened for 1.5s. Water drop occurs somewhere in there
    disp('Getting stimulus info from timestamp file');
    timestamps = readmatrix(timestamp_file);
    trials = unique(timestamps(:, 4));
    audio_tone{count} = zeros(size(fll_move{count}));
    water_drop{count} = zeros(size(fll_move{count}));
    for i = 1:length(trials)
        if trials(i) == 0, continue; end
        trial_start  = find(timestamps(:,4)==trials(i), 1); 
        audio_tone{count}(trial_start:trial_start+fs) = 1;
        if count > 2 
            % first 2 trials for each mouse are spontaneous trials - no
            % water drop
            water_drop{count}(trial_start+ceil(2*fs): trial_start+ceil(3.5*fs))=1;
        end
    end
    
    clear trials timestamps 
    count = count + 1;
%     if j == 3, break; end
%     master_SVD_file = [mouse_root_dir filesep 'masterSVD.mat'];
%     if ~isfile(master_SVD_file)
end