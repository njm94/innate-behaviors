% 

clear, close all
addpath('C:\Users\user\Documents\Nick\ridgeModel');
addpath('C:\Users\user\Documents\Nick\ridgeModel\widefield')
addpath('C:\Users\user\Documents\Nick\ridgeModel\smallStuff') 
addpath('C:\Users\user\Documents\Nick\grooming\utils')

%%
fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
for j = 1:7%length(data_list{1})
    data_dir = data_list{1}{j};
    disp(['Starting ' data_dir])
    fPath = [data_dir filesep 'ridge_outputs_video' filesep];
    if ~isdir(fPath), mkdir(fPath); end
    if isfile([fPath, 'summary.fig'])
        disp([data_dir, ' already processed. Skipping...'])
        continue
    end

    
    exp_date = strfind(data_dir, filesep);
    exp_date = data_dir(exp_date(end)+1:end);
    fs = 90;

    try
        brain_file = [data_dir, filesep, getAllFiles(data_dir, 'cam0_svd')];
        beh_file = [data_dir, filesep, getAllFiles(data_dir, [exp_date, '_svd'])];
        ME_file = [data_dir, filesep, get_file_with_str(data_dir, 'MEsvd')];
%         frame_file = [data_dir, filesep, get_file_with_str(data_dir, 'singleFrame')];
%         grooming_file = [data_dir, filesep, 'grooming_events_roi_filtered.mat'];
%         dlc_speed_file = [data_dir, filesep, getAllFiles(data_dir, 'speed.csv')];
        timestamp_file = [data_dir, filesep, getAllFiles(data_dir, 'trim.txt')];
%         snippets_dir = [data_dir filesep 'snippets'];
    catch
        disp('Missing one or more of the required files. Skipping...')
        continue
    end


    %% load brain svd components, multiply s into V
        % highpass filter
    disp('Loading brain data...')
    load(brain_file);
    Ubrain = U;

    % light-OFF signal may be captured in the brain data, resulting in a
    % massive filter artifact - remove a few frames just in case
    trial_length = size(V, 2) - 5; 
    Vbrain = s*V(:, 1:trial_length);

%     [b, a] = butter(2, 0.01/(fs/2), 'high');
    [b, a] = butter(1, [0.01 10]/(fs/2)); % bandpass
    Vbrain = filtfilt(b, a, Vbrain')';
    
    
%     %% load behavior svd components, multiply s into V
    load(beh_file)
    Ubeh = U;
    Vbeh = s*V(:, 1:trial_length);
    Vbeh = filtfilt(b, a, Vbeh')';
%     
%     %% load motion svd components, multiply s into V
%     load(ME_file)
%     Ume = U;
%     Vme = s*V(:, 1:trial_length);
    
    
    %%
    clear U s V
    

    %% load stimulus info from timestamp file
    % 10kHz tone occurs for 1 second immediately at start of trial
    disp('Getting audio stimulus info from timestampfile');
    timestamps = readmatrix(timestamp_file);
    trials = unique(timestamps(:, 4));
    audio_tone = zeros(trial_length, 1);
    for i = 1:length(trials)
        if trials(i) == 0, continue; end
        trial_start  = find(timestamps(:,4)==trials(i), 1); 
        audio_tone(trial_start:trial_start+fs) = 1;
    end
    
    %% build design matrix
    disp('Building design matrix')
    bopts.frameRate = fs;
    bopts.sPostTime=round(fs*1);
    bopts.framesPerTrial = trial_length; % nr. of frames per trial
    bopts.folds = 10; %nr of folds for cross-validation
    
    % make design matrix for audio stimulus
    [sMat, aEventIdx] = makeDesignMatrix(audio_tone, 2, bopts);
    moveLabels = {'Audio'};

    % concatenate the design matrix with temporal components from video and
    % motion energy
    fullR = [sMat Vbeh'];    
    
    % regressor index
    regIdx = [aEventIdx;  repmat(max(aEventIdx)+1, size(Vbeh,1), 1)];
    regLabels = [moveLabels, {'Video'}];
    
    
    %% run cross-validation
    %full model - this will take a moment
    disp('Running ridge regression with 10-fold cross-validation')
    [Vfull, fullBeta, ~, fullIdx, fullRidge, fullLabels] = crossValModel(fullR, Vbrain, regLabels, regIdx, regLabels, bopts.folds);
    save([fPath, char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_cvFull.mat'], 'Vfull', 'fullBeta', 'fullR', 'fullIdx', 'fullRidge', 'fullLabels', '-v7.3'); %save some results
    
    fullMat = modelCorr(Vbrain,Vfull,Ubrain) .^2; %compute explained variance
    
    %% run reduced models for unique contribution
    disp('Running reduced models')
    fig_flag = 1;
    reducedMat = [];
    for i = 1:length(regLabels)
        reduced = fullR;
        cIdx = regIdx == i;
%         if ~any(cIdx)
%             disp(['No ', regLabels{i}, '  events found. Skipping...'])
%             continue
%         end
        reduced(:, cIdx) = reduced(randperm(size(reduced, 1)), cIdx);
        
        [Vreduced{i}, reducedBeta{i}, reducedR, reducedIdx, reducedRidge, reducedLabels] = crossValModel(reduced, Vbrain, regLabels, regIdx, regLabels, bopts.folds);
        reducedMat(:, i) = modelCorr(Vbrain, Vreduced{i}, Ubrain) .^2; %compute explained variance
    end
    save([fPath, char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_cvReduced.mat'], 'Vreduced', 'reducedBeta', 'reducedR', 'reducedIdx', 'reducedRidge', 'reducedLabels', '-v7.3'); %save some results
    
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
            subplot(1, 3, i)
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

    savefig(gcf, [fPath 'summary.fig'])
end
  %%


