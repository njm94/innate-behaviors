
clear, close all, clc
% cd('/media/user/teamshare/nick/behavior/grooming/code/')
% addpath('/media/user/teamshare/nick/behavior/grooming/code/ridgeModel')
% addpath('/media/user/teamshare/nick/behavior/grooming/code/ridgeModel/widefield')
% addpath('/media/user/teamshare/nick/behavior/grooming/code/ridgeModel/smallStuff')

addpath('C:\Users\user\Documents\Nick\ridgeModel');
addpath('C:\Users\user\Documents\Nick\ridgeModel\widefield')
addpath('C:\Users\user\Documents\Nick\ridgeModel\smallStuff') 
addpath('C:\Users\user\Documents\Nick\grooming\utils');

%%

fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
groom_list = textscan(fileID, formatSpec);

fileID = fopen('expt3_datalist.txt','r');
formatSpec = '%s';
lick_list = textscan(fileID, formatSpec);

% data_list = [groom_list{1}; lick_list{1}];
data_list = groom_list{1};
current_mouse = '';


behaviors = {'Audio', 'Drop Left', 'Drop Right', 'Drop Center', 'Lick', ...
    'Elliptical', 'Elliptical Asymmetric', 'Left', 'Left Asymmetric', ...
    'Right', 'Right Asymmetric', 'Large Bilateral', 'FLL', 'FLR'};

grooming_behaviors = {'Elliptical', 'Elliptical Asymmetric', 'Left', ...
    'Left Asymmetric', 'Right', 'Right Asymmetric', 'Large Bilateral'};

% expts_to_analyze = ;
thy1_idx = 1:7;
ai94_idx = 8:13;
camk_idx = 14:25;
%%

for j = ai94_idx(1):ai94_idx(end)+1%:20%_file1:length(data_list)+1
     try
        data_dir = data_list{j};
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

                % regressors
                audio_tone{i} = audio_tone{i}(1:min_trial_length);
                drop_left{i} = drop_left{i}(1:min_trial_length);
                drop_right{i} = drop_right{i}(1:min_trial_length);

                lick_regressor{i} = lick_regressor{i}(1:min_trial_length);
                grooming_regressors{i} = grooming_regressors{i}(1:min_trial_length,:);

%                 left{i} = left{i}(1:min_trial_length);
%                 right{i} = right{i}(1:min_trial_length);
%                 elliptical{i} = elliptical{i}(1:min_trial_length);
%                 large_left{i} = largeleft{i}(1:min_trial_length);
%                 large_right{i} = largeright{i}(1:min_trial_length);
%                 large_bilateral{i} = largebilateral{i}(1:min_trial_length);
%                 lick_regressor{i} = lick_regressor{i}(1:min_trial_length);

                fll_move{i} = fll_move{i}(1:min_trial_length);
                flr_move{i} = flr_move{i}(1:min_trial_length);

            end
            Vmaster = catcell(2, Vmaster);
            audio_tone = catcell(1, audio_tone);
            drop_left = catcell(1, drop_left);
            drop_right = catcell(1, drop_right);


            audio_tone = get_on_time(audio_tone);
            drop_left = get_on_time(drop_left);
            drop_right = get_on_time(drop_right);

%             left = catcell(2, left)';
%             right = catcell(2, right)';
%             elliptical = catcell(2, elliptical)';
%             large_left = catcell(2, large_left)';
%             large_right = catcell(2, large_right)';
%             large_bilateral = catcell(2, large_bilateral)';
            grooming_regressors = catcell(1, grooming_regressors);
            lick_regressor = catcell(1, lick_regressor);
%             [lick_start, lick_timer] = start_timer(lick_regressor, 3, fs); 

%             % Aggregate regressors and create timers
%             bilateral_grooming = any(grooming_regressors(:, [1,4,6,7]), 2);
%             [bilat_start_idx, bilat_timer] = start_timer(bilateral_grooming, 0.25, fs); 
% 
%             all_groom = any(grooming_regressors,2);
%             [groom_start, groom_timer] = start_timer(all_groom, 3, fs); 

            fll_move = catcell(1, fll_move);
            flr_move = catcell(1, flr_move);
%%
            [~, activity_timer] = start_timer(any([grooming_regressors, lick_regressor, fll_move, flr_move],2), 2, fs);
%%
%             fll_start = start_timer(fll_move, 3, fs);
%             flr_start = start_timer(flr_move, 3, fs);

            % ridge here
            disp('Building design matrix')
            opts.frameRate = fs;
            opts.sPostTime=round(fs*1);
            opts.mPreTime = ceil(0.5 * fs);  % precede motor events to capture preparatory activity in frames (used for eventType 3)
            opts.mPostTime = ceil(2 * fs);   % follow motor events for mPostStim in frames (used for eventType 3)
            opts.framesPerTrial = min_trial_length; % nr. of frames per trial
            opts.folds = 10; %nr of folds for cross-validation
            
%             regressor_mat = [audio_tone' water_drop' ipsi' contra' bilat' fll_move' flr_move']; % tone drop ipsi contra bilatInclusive forelimbInc
            regressor_mat = [audio_tone drop_left drop_right grooming_regressors lick_regressor fll_move flr_move];
%             regressor_mat = [audio_tone drop_left drop_right elliptical large_left large_right large_bilateral left right lick_regressor fll_move flr_move]; % tone drop ipsi contra bilatInclusive forelimbInc

            % Full-Trial events:    new_trial
            % Post-Stimulus events: audio_tone, water_drop
            % Peri-Stimulus events: ipsi, contra, bilat, fll_move, flr_move
            [dMat, regIdx] = makeDesignMatrix(regressor_mat, [2, 2, 2, repmat(3, 1, size(grooming_regressors,2)), 3, 3, 3], opts);
%             [dMat, regIdx] = makeDesignMatrix(regressor_mat, [2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3], opts);
%             regLabels = {'audio', 'drop', 'ipsi', 'contra', 'bilat', 'left_move', 'right_move'}; %some movement variables
%             regLabels = {'Audio', 'DropLeft', 'DropRight', ...
%                 'Elliptical', 'LargeLeft', 'LargeRight', 'LargeBilateral', ...
%                 'Left', 'Right', ...
%                 'Lick', 'LeftMove', 'RightMove'};

            regLabels = [{'Audio', 'DropLeft', 'DropRight'}, ...
                grooming_behaviors, ...
                {'Lick', 'LeftMove', 'RightMove', 'Timer'}];

            regIdx = [regIdx; max(regIdx)+1;];
            
            activity_timer = activity_timer ./ fs; % scale the activity timer so it's expressed wrt seconds
            fullR = [dMat, activity_timer];

            disp('Running ridge regression with 10-fold cross-validation')
            [Vfull, fullBeta, ~, fullIdx, fullRidge, fullLabels] = crossValModel(fullR, Vmaster, regLabels, regIdx, regLabels, opts.folds);
            save([fPath, char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_cvFull.mat'], 'Vfull', 'fullBeta', 'fullR', 'fullIdx', 'fullRidge', 'fullLabels', '-v7.3'); %save some results
            
            fullMat = modelCorr(Vmaster,Vfull,Umaster) .^2;

            break

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

            save([fPath, char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_cvReduced.mat'], 'Vreduced', 'reducedBeta', 'reducedR', 'reducedIdx', 'reducedRidge', 'reducedLabels', '-v7.3'); %save some results

            figure
            subplot(5, 3, 1)
            imagesc(reshape(fullMat, [128 128]))
            xticks([])
            c=colorbar;
            title('FullMat')
            c.Label.String = 'cvR^2';
            yticks([])
            for i = 1:length(regLabels)
                subplot(5, 3, i+1)
                imagesc(reshape(fullMat - reducedMat(:,i), [128 128]))
                c=colorbar;
                title(regLabels{i})
                c.Label.String = '\DeltaR^2';
            
                xticks([])
                yticks([])
            end

        
            savefig(gcf, [fPath, char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_summary.fig'])
            break

            % all mice completed - break the loop
            if j == length(data_list{1})+1, break; end
        end        
        count = 1;
        disp('Loading master basis set')
        load([mouse_root_dir filesep 'Umaster.mat'])
        clear Vmaster left right elliptical large_left large_right largebilateral lick_regressor fll_move flr_move audio_tone drop_left drop_right
        current_mouse = mouse_id;

        fPath = [mouse_root_dir filesep 'outputs' filesep];

        %     %% draw mask - might be needed for memory management
        %     frame = loadtiff(frame_file);
        %     mask = draw_roi(frame, 4);
        
    end

    brain_file = [data_dir, filesep, getAllFiles(data_dir, 'cam0_svd')];
    dlc_speed_file = [data_dir, filesep, getAllFiles(data_dir, 'speed.csv')];
    timestamp_file = [data_dir, filesep, getAllFiles(data_dir, 'trim.txt')];
    boris_file = [data_dir, filesep, getAllFiles(data_dir, 'events.tsv')];
    config_file = [data_dir, filesep, getAllFiles(data_dir, 'config.ini')];

    if ~any(isfile(brain_file) | isfile(dlc_speed_file) | isfile(timestamp_file) | isfile(boris_file))
        disp('Missing one or more of the required files. Skipping...')
        continue        
    end

    % read ini file
    ini = ini2struct(config_file);
    fpath = fieldnames(ini);
    ini = ini.(fpath{contains(fpath, 'sentech')});
    fs = str2double(ini.framerate);
    
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
    Vbrain = recastV(Umaster, U, s, V(:, 1:trial_length));
    % Vmaster{count} = Vbrain;
    [b, a] = butter(2, 0.01/(fs/2), 'high');
    Vmaster{count} = filtfilt(b, a, Vbrain')';
    clear U s V Vbrain

    %% load grooming events
    disp('Reading BORIS')
    [events, b_idx, b_tab] = read_boris(boris_file);
    video_end = find(events.('Video End'));
    vid_end_idx = contains(events.Properties.VariableNames, 'Video End');
    events = removevars(events, vid_end_idx);
    b_idx(vid_end_idx) = [];

    % resample brain to match video_end size
    Vmaster{count} = resamplee(Vmaster{count}', video_end, trial_length)';

    % consolidate lick events and remove from table
    lick_idx = find(contains(events.Properties.VariableNames, 'Lick'));
    lick_regressor{count} = zeros(video_end, 1);
    for i = 1:length(lick_idx)
        lick_regressor{count}(b_idx{lick_idx(i)}(:,1)) = 1;
    end
    events = removevars(events, lick_idx);
    b_idx(lick_idx) = [];

    % consolidate drop events
    drop_idx = contains(events.Properties.VariableNames, 'Drop');
    if any(drop_idx)
        events = removevars(events, drop_idx);
        b_idx(drop_idx) = [];
    end

    % remaining events are grooming events
    grooming_regressors{count} = zeros(video_end, length(grooming_behaviors));
    for i = 1:length(grooming_behaviors)
        if any(strcmpi(events.Properties.VariableNames, grooming_behaviors{i}))
            idx = find(strcmpi(events.Properties.VariableNames, grooming_behaviors{i}));
            grooming_regressors{count}(b_idx{idx}(:,1),i) = 1;
        end
    end


    %% load DLC tracks
    disp('Loading DLC tracks')
    dlc_speed = readmatrix(dlc_speed_file);
    fll_speed = dlc_speed(1:video_end,1);
    flr_speed = dlc_speed(1:video_end,2);
    
    % binarize movement speed
    fll_speed = fll_speed > mean(fll_speed) + std(fll_speed);
    flr_speed = flr_speed >  mean(flr_speed) + std(flr_speed);
    
%     k = 35; % smoothing kernel - this matches smoothing kernel in grooming detection but may need to change
    fll_move{count} = get_on_time(fll_speed);
    flr_move{count} = get_on_time(flr_speed);


    clear dlc_speed fll_speed flr_speed

    %% load stimulus info from timestamp file
    % 10kHz tone occurs for 1 second immediately at start of trial
    % Tone is followed by 1 second of silence
    % Then valve is opened for 1.5s. Water drop occurs somewhere in there
    disp('Getting stimulus info from timestamp file');
    timestamps = readmatrix(timestamp_file);
    trials = unique(timestamps(:, 4));
    audio_tone{count} = zeros(size(fll_move{count}));
    drop_left{count} = zeros(size(fll_move{count}));
    drop_right{count} = zeros(size(fll_move{count}));
    for i = 1:length(trials)
        if trials(i) == 0, continue; end
        trial_start  = find(timestamps(:,4)==trials(i), 1); 
        audio_tone{count}(trial_start:trial_start+fs) = 1;
        if count > 2 && count < 6
            % first 2 trials for each mouse are spontaneous trials - no
            % water drop. then on last trial for HYL3, drop goes on left
            % side
            drop_right{count}(trial_start+ceil(2*fs): trial_start+ceil(3.5*fs))=1;
        elseif count == 6
            drop_left{count}(trial_start+ceil(2*fs): trial_start+ceil(3.5*fs))=1;
        end
    end
    
    clear trials timestamps 
    count = count + 1;
%     if j == 3, break; end
%     master_SVD_file = [mouse_root_dir filesep 'masterSVD.mat'];
%     if ~isfile(master_SVD_file)



% %% load behavior video SVD

end


function [start_idx, behavior_timer] = start_timer(data, k, fs)
% takes a binary event vector (data) and aggregation window size (k) in
% seconds
% aggregates all behavior events and creates a linear ramp
% also outputs start of each aggregated behavior session


idx = arr2idx(aggregate(data, k, fs));
start_idx = zeros(size(data));
start_idx(idx(:,1)) = 1;

behavior_timer = zeros(size(data));
for i = 1:size(idx,1)
    tmp = 1:(idx(i,2)-idx(i,1))+1;
    behavior_timer(idx(i,1):idx(i,2)) = tmp;
end

end


function new_dat = get_on_time(data)
new_dat = zeros(size(data));
d = diff(data);
t_on = find(d>0)+ 1;

new_dat(t_on) = 1;

end