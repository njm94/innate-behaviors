% 

clear,
% close all
addpath('C:\Users\user\Documents\Nick\ridgeModel');
addpath('C:\Users\user\Documents\Nick\ridgeModel\widefield')
addpath('C:\Users\user\Documents\Nick\ridgeModel\smallStuff') 
addpath('C:\Users\user\Documents\Nick\grooming\utils')

%%
fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
data_list = textscan(fileID, formatSpec);

count = 1;
for j = 7%1:length(data_list{1})

    for num_comps = 10:10:1000

    data_dir = data_list{1}{j};
    disp(['Starting ' data_dir])
    fPath = [data_dir filesep 'ridge_outputs_video' filesep];
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

        % THIS IS A CONTROL TO SEE IF WE CAN REPRODUCE THE BRAIN DATA AFTER
        % SHUFFLING THE BEHAVIOR VIDEOS
%         ME_file = 'Y:\nick\behavior\grooming\1p\ECR2_thy1\20231116161015\ECR2_thy1_20231116161015_MEsvd.mat';
%         beh_file = 'Y:\nick\behavior\grooming\1p\IBL2_tTA\20231120033323\IBL2_tTA_20231120033323_svd.mat';
        % -----------------------------------

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
    Ubrain = U(:,1:num_comps);

    % light-OFF signal may be captured in the brain data, resulting in a
    % massive filter artifact - remove a few frames just in case
    trial_length = size(V, 2) - 5; 
    Vbrain = s*V(:, 1:trial_length);
    Vbrain = Vbrain(1:num_comps,:);

%     [b, a] = butter(2, 0.01/(fs/2), 'high');
    [b, a] = butter(1, [0.01 10]/(fs/2)); % bandpass
    Vbrain = filtfilt(b, a, Vbrain')';
    
    
%     %% load behavior svd components, multiply s into V
    load(beh_file)
    Ume = U;
    Vme = s*V(:, 1:trial_length);
%     Vme = filtfilt(b, a, Vme')';
%     
%     %% load motion svd components, multiply s into V
%     load(ME_file)
%     Ume = U;
%     Vme = s*V(:, 1:trial_length);
%     Vme = filtfilt(b, a, Vme')';
    
    
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
%         audio_tone(trial_start:trial_start+fs) = 1;
        audio_tone(trial_start) = 1;
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
    
%     fullR = [sMat Vbeh'];   
    fullR = [sMat Vme(1:num_comps,:)'];   
%     fullR = Vme(1:num_comps,:)';
    
    % regressor index
    regIdx = [aEventIdx;  repmat(max(aEventIdx)+1, num_comps, 1)];
    regLabels = [moveLabels, {'MotionEnergy'}];

% %% run QR and check for rank-defficiency. This will show whether a given regressor is highly collinear with other regressors in the design matrix.
% % The resulting plot ranges from 0 to 1 for each regressor, with 1 being
% % fully orthogonal to all preceeding regressors in the matrix and 0 being
% % fully redundant. Having fully redundant regressors in the matrix will
% % break the model, so in this example those regressors are removed. In
% % practice, you should understand where the redundancy is coming from and
% % change your model design to avoid it in the first place!
% 
% rejIdx = false(1,size(fullR,2));
% [~, fullQRR] = qr(bsxfun(@rdivide,fullR,sqrt(sum(fullR.^2))),0); %orthogonalize normalized design matrix
% figure; plot(abs(diag(fullQRR)),'linewidth',2); ylim([0 1.1]); title('Regressor orthogonality'); drawnow; %this shows how orthogonal individual regressors are to the rest of the matrix
% axis square; ylabel('Norm. vector angle'); xlabel('Regressors');
% if sum(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1))) < size(fullR,2) %check if design matrix is full rank
%     temp = ~(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1)));
%     fprintf('Design matrix is rank-defficient. Removing %d/%d additional regressors.\n', sum(temp), sum(~rejIdx));
%     rejIdx(~rejIdx) = temp; %reject regressors that cause rank-defficint matrix
% end
% % save([fPath filesep 'regData.mat'], 'fullR', 'regIdx', 'regLabels','fullQRR','-v7.3'); %save some model variables

    
    %% run cross-validation
    %full model - this will take a moment
    disp('Running ridge regression with 10-fold cross-validation')
    [Vfull, fullBeta, ~, fullIdx, fullRidge, fullLabels] = crossValModel(fullR, Vbrain, regLabels, regIdx, regLabels, bopts.folds);
%     save([fPath, char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), 'movie_cvFull_', num2str(num_comps),'.mat'], 'Vfull', 'fullBeta', 'fullR', 'fullIdx', 'fullRidge', 'fullLabels', '-v7.3'); %save some results
    
    fullMat(:,:,count) = modelCorr(Vbrain,Vfull,Ubrain(:,1:num_comps)) .^2; %compute explained variance
    count = count + 1;
    end
%     %% run reduced models for unique contribution
%     disp('Running reduced models')
%     fig_flag = 1;
%     reducedMat = [];
%     for i = 1:length(regLabels)
%         reduced = fullR;
%         cIdx = regIdx == i;
% %         if ~any(cIdx)
% %             disp(['No ', regLabels{i}, '  events found. Skipping...'])
% %             continue
% %         end
%         reduced(:, cIdx) = reduced(randperm(size(reduced, 1)), cIdx);
%         
%         [Vreduced{i}, reducedBeta{i}, reducedR, reducedIdx, reducedRidge, reducedLabels] = crossValModel(reduced, Vbrain, regLabels, regIdx, regLabels, bopts.folds);
%         reducedMat(:, i) = modelCorr(Vbrain, Vreduced{i}, Ubrain) .^2; %compute explained variance
%     end
%     save([fPath, char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), 'movie_cvReduced.mat'], 'Vreduced', 'reducedBeta', 'reducedR', 'reducedIdx', 'reducedRidge', 'reducedLabels', '-v7.3'); %save some results
    
    %%
    
%     figure, 
%     if fig_flag
%         subplot(1, length(regLabels)+1, 1)
%         imagesc(reshape(fullMat, [128 128])')
%         xticks([])
%         c=colorbar;
%         title('FullMat')
%         c.Label.String = 'cvR^2';
%         yticks([])
%         for i = 1:length(regLabels)
%             subplot(1, 3, i)
%             imagesc(reshape(fullMat - reducedMat(:,i), [128 128])')
%             c=colorbar;
%             title(regLabels{i})
%             c.Label.String = '\DeltaR^2';
%         
%             xticks([])
%             yticks([])
%         end
%     else
%         imagesc(reshape(fullMat, [128 128])')
%         xticks([])
%         c=colorbar;
%         title('FullMat')
%         c.Label.String = 'cvR^2';
%         yticks([])
%     end
% 
%     savefig(gcf, [fPath 'movie_summary.fig'])
end
  %%

testreal = permute(reshape(Ubrain(:,1:num_comps) * Vbrain(:,6.45e4:7.1e4), 128, 128, []), [2 1 3]);
testmodel = permute(reshape(Ubrain(:,1:num_comps) * Vfull(:,6.45e4:7.1e4), 128, 128, []), [2 1 3]);

resid = testreal-testmodel;


%%

bvid = permute(reshape(Ume(:,1:num_comps)*Vme(1:num_comps, 6.45e4:7.1e4), 214, 107, []), [2 1 3]);
bvidall = permute(reshape(Ume*Vme(:, 6.45e4:7.1e4), 214, 107, []), [2 1 3]);


%%

figure
show_mov(zscore([padarray(bvid, 21, 0, 'pre'), resid], [], 3))

%%
N = size(bvid,3);
    for i = 1:N
    ff = figure(1);
    ff.Position= [99 502 821 311];
    subplot(1,2,1)
    imagesc(bvid(:,:,i))
    axis off
    clim([0 100])
    colormap gray
    freezeColors()

    subplot(1,2,2)
    imagesc(resid(:,:,i))
    axis off
    colorbar
    
    clim([-200 200])
    colormap(bluewhitered())
    freezeColors()

      F(i) = getframe(gcf) ;
      drawnow
    end
  % create the video writer with 1 fps
  writerObj = VideoWriter('Y:\nick\behavior\grooming\1p\IBL2_tTA\20231120033323\myVideo.avi');
  writerObj.FrameRate = 60;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

  %% visualize the reconstruction




fullLabels
visual = true;
cBetaRight = check_beta('Audio', fullLabels, fullIdx, Ubrain, fullBeta{1}, Vfull, [], visual);

%% visulize pixel weights

test = fullBeta{1}(92:end,:);


UU = Ume(:,1:200) * test;


