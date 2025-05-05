% 

clear, close all

if ~isunix
    addpath('C:\Users\user\Documents\Nick\ridgeModel');
    addpath('C:\Users\user\Documents\Nick\ridgeModel\widefield')
    addpath('C:\Users\user\Documents\Nick\ridgeModel\smallStuff') 
    addpath('C:\Users\user\Documents\Nick\grooming\utils')
else
    addpath('/home/user/Documents/grooming/ridgeModel');
    addpath('/home/user/Documents/grooming/ridgeModel/widefield')
    addpath('/home/user/Documents/grooming/ridgeModel/smallStuff') 
    addpath('/home/user/Documents/grooming/utils')
end


%%
fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
for j = 7%length(data_list{1})
    data_dir = data_list{1}{j};
       
    if isunix
        data_dir = strrep(data_dir, '\', '/');
        data_dir = strrep(data_dir, 'Y:', '/media/user/teamshare');
    end
    disp(['Starting ' data_dir])

    fPath = [data_dir filesep 'ridge_outputs' filesep];
    if ~isdir(fPath), mkdir(fPath); end
%     if isfile([fPath, 'summary.fig'])
%         disp([data_dir, ' already processed. Skipping...'])
%         continue
%     end

    
    exp_date = strfind(data_dir, filesep);
    exp_date = data_dir(exp_date(end)+1:end);
    fs = 90;

    try
        brain_file = [data_dir, filesep, getAllFiles(data_dir, 'cam0_svd')];
        beh_file = [data_dir, filesep, getAllFiles(data_dir, [exp_date, '_svd'])];
        ME_file = [data_dir, filesep, get_file_with_str(data_dir, 'MEsvd')];
%         frame_file = [data_dir, filesep, get_file_with_str(data_dir, 'singleFrame')];
%         grooming_file = [data_dir, filesep, 'grooming_events_roi_filtered.mat'];
        dlc_speed_file = [data_dir, filesep, getAllFiles(data_dir, 'speed.csv')];
        timestamp_file = [data_dir, filesep, getAllFiles(data_dir, 'trim.txt')];
%         snippets_dir = [data_dir filesep 'snippets'];
        boris_file = [data_dir, filesep, get_file_with_str(data_dir, 'events.tsv')];
    catch
        disp('Missing one or more of the required files. Skipping...')
        continue
    end

%     if any([~isfile(brain_file), ~isfile(dlc_speed_file), ~isfile(beh_file), ~isfile(ME_file), ~isfile(grooming_file)])
%         disp('Missing one or more of the required files. Skipping...')
%         continue
%     end
    
    %% load grooming events from BORIS file
    [events, b_idx, t] = read_boris(boris_file);

    %% load brain svd components, multiply s into V
        % highpass filter
    disp('Loading brain data...')
    load(brain_file);
    Ubrain = U;

    [b, a] = butter(1, [0.01 10]/(fs/2)); % bandpass

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
%     Vbrain = recastV(Umaster, U, s, V(:, 1:trial_length));
    Vbrain = s*V;
    Vbrain = filtfilt(b, a, Vbrain(:,1:trial_length)')';

%     % light-OFF signal may be captured in the brain data, resulting in a
%     % massive filter artifact - remove a few frames just in case
%     trial_length = size(V, 2) - 5; 
%     Vbrain = s*V(:, 1:trial_length);

%     [b, a] = butter(2, 0.01/(fs/2), 'high');

%     Vbrain = filtfilt(b, a, Vbrain')';
    
    
% %     %% load behavior svd components, multiply s into V
%     load(beh_file)
%     Ubeh = U;
%     Vbeh = s*V;
% %     
% %     %% load motion svd components, multiply s into V
%     load(ME_file)
%     Ume = U;
%     Vme = s*V;
    
    
    %%
    clear U s V
    

    %% 
    [behaviors, labels] = parse_snippets(snippets_dir);
    bmat = zeros(length(behaviors), size(Vbrain,2));
    for i = 1:length(behaviors)
        for j = 1:size(behaviors{i}, 1)
            bmat(i, behaviors{i}(j,1):behaviors{i}(j,2)) = 1;
        end
    end

    %% load DLC tracks
    disp('Loading DLC tracks')
    dlc_speed = readmatrix(dlc_speed_file);
    fll_speed = dlc_speed(:,1);
    flr_speed = dlc_speed(:,2);

    LeftMove = fll_speed > mean(fll_speed) + std(fll_speed);
    RightMove = flr_speed >  mean(flr_speed) + std(flr_speed);

    LeftMove = movmax(LeftMove, 35);
    RightMove = movmax(RightMove, 35);


    %% load stimulus info from timestamp file
    % 10kHz tone occurs for 1 second immediately at start of trial
    % Tone is followed by 1 second of silence
    % Then valve is opened for 1.5s. Water drop occurs somewhere in there
    disp('Getting stimulus info from timestampfile');
    timestamps = readmatrix(timestamp_file);
    trials = unique(timestamps(:, 4));
    Audio = zeros(size(LeftMove));
%     water_drop = zeros(size(LeftMove));
    for i = 1:length(trials)
        if trials(i) == 0, continue; end
        trial_start  = find(timestamps(:,4)==trials(i), 1); 
        Audio(trial_start:trial_start+fs) = 1;
%         water_drop(trial_start+ceil(2*fs): trial_start+ceil(3.5*fs))=1;
    end

    bvars = ["Drop Hits Left", "Drop Hits Right", "Drop Hits Center", "Audio", ...
        "Lick Left", "Lick Right", "Lick Unknown", "Lick No Paws", ...
        "Right", "Left", "Right Asymmetric", "Left Asymmetric" ...
        "Elliptical", "Elliptical Asymmetric", "Large Bilateral", ...
        "LeftMove", "RightMove"];


    
    %% build design matrix
    disp('Building design matrix')
    bopts.frameRate = fs;
    bopts.sPostTime=round(fs*1);
    bopts.mPreTime = ceil(0.5 * fs);  % precede motor events to capture preparatory activity in frames (used for eventType 3)
    bopts.mPostTime = ceil(2 * fs);   % follow motor events for mPostStim in frames (used for eventType 3)
    bopts.framesPerTrial = trial_length; % nr. of frames per trial
    bopts.folds = 10; %nr of folds for cross-validation

    % consolidate lick events
    lick_idx = find(contains(events.Properties.VariableNames, 'Lick'));
    lick_events = zeros(trial_length, 1);
    for pp = 1:length(lick_idx)
        lick_events(b_idx{lick_idx(pp)}(:,1)) = 1;
    end


    
    % remove lick and point events from event matrix
    idx = contains(events.Properties.VariableNames, 'Lick') | ...
        contains(events.Properties.VariableNames, 'Video');
    stroke_events = removevars(events, idx);


    b_idx2 = b_idx(~idx);


%     bmat = table2array(events(:,~strcmp(events.Properties.VariableNames, 'Video End')));
%     bmat = [Audio(1:trial_length) bmat(1:trial_length, :)];

    bmat = zeros(size(stroke_events));

    % specify the type of variable for making the design matrix
    reg_type = [];
    for ii = 1:size(stroke_events,2)
        if contains(events.Properties.VariableNames{ii}, 'Drop')
            reg_type = [reg_type 2];
            bmat(:,ii) = table2array(stroke_events(:,ii));
        else
            reg_type = [reg_type 3];
            bmat(b_idx2{ii}(:,1), ii) = 1;
        end
    end

    % concatenate motor variable type for lick
    bmat = [bmat(1:trial_length,:), lick_events];
    reg_type = [reg_type, 3];


    % concatenate another STIMULUS variable type for the Audio tone
    bmat = [Audio(1:trial_length) bmat(1:trial_length, :)];
    reg_type = [2 reg_type];

    

%     bmat = [FL_R' FL_L']; %ipsi contra
%     bmat = [ipsi' contra' bilat']; % ipsi contra bilatExclusive
%     bmat = [FL_R' FL_L' bilat']; % ipsi contra bilatIcnlusive
%     bmat = [ipsi' contra' fll_movement flr_movement];
%     bmat = [FL_R' FL_L' bilat', fll_movement_ex, flr_movement_ex]; % ipsi contra bilatInclusive forelimbEx
%     bmat = [audio_tone water_drop FL_R' FL_L' lick2' fll_movement flr_movement]; % tone drop ipsi contra bilatInclusive forelimbInc
%    bmat = [Audio water_drop bmat', LeftMove, RightMove];

%     [movMat, movEventIdx2] = makeDesignMatrix(bmat, [3, 3, 3, 3, 3], bopts);
%     [movMat, movEventIdx2] = makeDesignMatrix(bmat, [3, 3], bopts);
%     [movMat, movEventIdx2] = makeDesignMatrix(bmat, [3, 3, 3], bopts);
    [movMat, regIdx] = makeDesignMatrix(bmat, reg_type, bopts);

%     moveLabels = {'audio', 'drop', 'ipsi', 'contra', 'lick' 'left_move', 'right_move'}; %some movement variables
%     moveLabels = {'ipsi', 'contra', 'left_move', 'right_move'}; %some movement variables    
%     moveLabels = {'ipsi', 'contra', 'bilat', 'left_move', 'right_move'}; %some movement variables
%     moveLabels = {'ipsi', 'contra'}; %some movement variables
%     moveLabels = {'ipsi', 'contra', 'bilat'}; %some movement variables
% moveLabels = {'audio', 'drop', 'lick', 'elliptical', 'largeleft', 'largeright', 'largebilateral', 'right', 'left', 'leftmove', 'rightmove'};

% moveLabels = {'Audio', 'Drop', 'Lick', ...
%                 'Elliptical', 'LargeLeft', 'LargeRight', 'LargeBilateral', ...
%                 'Right', 'Left', ...
%                 'LeftMove', 'RightMove'};

regLabels = stroke_events.Properties.VariableNames;
regLabels = ['Audio', regLabels, 'Lick'];

    % fullR = [Vbeh' Vme'];
    % fullR = [Vme'];    
    fullR = [movMat fll_speed(1:trial_length), flr_speed(1:trial_length)];
    regIdx = [regIdx; max(regIdx)+1; max(regIdx)+2];
    regLabels = [regLabels, 'FLL', 'FLR'];
%     fullR = [fullR Vbeh(:,1:size(fullR,1))'];    
    
%     %% run ridge
%     [ridgeVals, dimBeta] = ridgeMML(Vbrain', fullR, true); %get ridge penalties and beta weights.
%     Vmodel = (fullR * dimBeta)';
%     %%
%     corrMat = reshape(modelCorr(Vbrain,Vmodel,Ubrain) .^2, [128 128])'; %compute explained variance
%     figure, 
%     imagesc(mask .* corrMat), 
%     c=colorbar;
%     c.Label.String = 'cvR^2';
%     c.FontSize =12;
    
    % % regIdx = [taskIdx; moveIdx + max(taskIdx); repmat(max(moveIdx)+max(taskIdx)+1, size(vidR,2), 1)]; %regressor index
%     regIdx = [movEventIdx2;  repmat(max(movEventIdx2)+1, size(Vbeh',2), 1)];
%     regLabels = [moveLabels, {'video'}];
    % 
    %%
    cVar = 'Lick';
    
    % find beta weights for current variable
    cIdx = regIdx == find(ismember(regLabels,cVar));
    U = reshape(Ubrain, [], size(Vbrain,1)); 
    cBeta = U * fullBeta{1}(cIdx, :)';
    cBeta = permute(reshape(cBeta, 128, 128, []), [2, 1, 3]); 
    U = reshape(U, 128, 128, size(Vbrain,1)); 
    show_mov(cBeta)
    % compareMovie(cBeta.*mask)
    
    
    %% run cross-validation
    
    bopts.folds = 1;

    %full model - this will take a moment
    disp('Running ridge regression with 10-fold cross-validation')
    [Vfull, fullBeta, ~, fullIdx, fullRidge, fullLabels] = crossValModel(fullR, Vbrain, regLabels, regIdx, regLabels, bopts.folds);
    save([fPath, char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_cvFull.mat'], 'Vfull', 'fullBeta', 'fullR', 'fullIdx', 'fullRidge', 'fullLabels', '-v7.3'); %save some results
    
    fullMat = modelCorr(Vbrain,Vfull,Ubrain) .^2; %compute explained variance
    
    %% run reduced models for unique contribution
%     if size(event_type,1) <= 1
%         fig_flag = 0;
%         disp('Not enough unique event types to run reduce models. Skipping...')
%     else
        disp('Running reduced models')
        fig_flag = 1;
        reducedMat = [];
        for i = 1:length(regLabels)
            reduced = fullR;
            cIdx = regIdx == i;
            if ~any(cIdx)
                disp(['No ', regLabels{i}, '  events found. Skipping...'])
                continue
            end
            reduced(:, cIdx) = reduced(randperm(size(reduced, 1)), cIdx);
        
            [Vreduced{i}, reducedBeta{i}, reducedR, reducedIdx, reducedRidge, reducedLabels] = crossValModel(reduced, Vbrain, regLabels, regIdx, regLabels, bopts.folds);
            reducedMat(:, i) = modelCorr(Vbrain, Vreduced{i}, Ubrain) .^2; %compute explained variance
        end
    %     reducedMat = reshape(reducedMat, 128, 128, length(regLabels));
        save([fPath, char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_cvReduced.mat'], 'Vreduced', 'reducedBeta', 'reducedR', 'reducedIdx', 'reducedRidge', 'reducedLabels', '-v7.3'); %save some results
%     end


%     %% run reduced models for unique contribution
%     if size(event_type,1) <= 1
%         fig_flag = 0;
%         disp('Not enough unique event types to run reduce models. Skipping...')
%     else
%         disp('Running reduced models')
%         fig_flag = 1;
%         reducedMat = [];
%         for i = 1:length(regLabels)
%             reduced = fullR;
%             cIdx = regIdx == i;
%             if ~any(cIdx)
%                 disp(['No ', regLabels{i}, '  events found. Skipping...'])
%                 continue
%             end
%             reduced(:, cIdx) = reduced(randperm(size(reduced, 1)), cIdx);
%         
%             [Vreduced{i}, reducedBeta{i}, reducedR, reducedIdx, reducedRidge, reducedLabels] = crossValModel(reduced, Vbrain, regLabels, regIdx, regLabels, bopts.folds);
%             reducedMat(:, i) = modelCorr(Vbrain, Vreduced{i}, Ubrain) .^2; %compute explained variance
%         end
%     %     reducedMat = reshape(reducedMat, 128, 128, length(regLabels));
%         save([fPath, char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_cvReduced.mat'], 'Vreduced', 'reducedBeta', 'reducedR', 'reducedIdx', 'reducedRidge', 'reducedLabels', '-v7.3'); %save some results
%     end

    
    %%
    
    figure, 
    if fig_flag
        subplot(1, length(regLabels)+1, 1)
        imagesc(reshape(fullMat, [128 128])')
        xticks([])
        c=colorbar;
        title('FullMat')
        c.Label.String = 'cvR^2';
        yticks([])
        for i = 1:length(regLabels)
            subplot(4, 3, i)
            imagesc(reshape(fullMat - reducedMat(:,i), [128 128])')
            c=colorbar;
            title(regLabels{i})
            c.Label.String = '\DeltaR^2';
        
            xticks([])
            yticks([])
        end
    else
        imagesc(reshape(fullMat, [128 128])')
        xticks([])
        c=colorbar;
        title('FullMat')
        c.Label.String = 'cvR^2';
        yticks([])
    end
%%
    savefig(gcf, [fPath 'summary.fig'])
    %%
end
  %%


