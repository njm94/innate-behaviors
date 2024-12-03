% Regress movement out of 2p data

clear, clc

% addpath(genpath('C:\Users\user\Documents\Nick\grooming'));
addpath('C:\Users\user\Documents\Nick\grooming\utils')
[event_file, event_path] = uigetfile('*.tsv','Select BORIS event labels.', 'Y:\nick\behavior\grooming\2p');
[events, b_idx, ~] = read_boris([event_path, filesep, event_file]);

if ~isfile([event_path, 'Nresample.mat'])
    [neuron_file, neuron_path] = uigetfile('*clean.mat','Select cleaned neuron data.', event_path);
    load([neuron_path, filesep, neuron_file]);
    
    disp("Resampling neural data to match behavior")
    Nresample = resample(N, size(events,1), size(N, 2), 'Dimension', 2);

    save([neuron_path, 'Nresample.mat'], 'Nresample', 'cstat', 'nloc', 'tforms')
else
    disp('Loading resampled neuron data')
    load([event_path, 'Nresample.mat'])
end


%% load DLC tracks
vel = readmatrix([event_path, filesep, getAllFiles(event_path, '_vel.csv')]);
vid_end = find(events.("Video End"));

vel = vel(1:vid_end,:);
flrv = sum(vel(:,4:5).^2, 2).^0.5;
fllv = sum(vel(:,7:8).^2, 2).^0.5;


%% build design matrix
fs = 90;
disp('Building design matrix')
bopts.frameRate = 90;
bopts.sPostTime=round(fs*3);
bopts.mPreTime = ceil(0.5 * fs);  % precede motor events to capture preparatory activity in frames (used for eventType 3)
bopts.mPostTime = ceil(2 * fs);   % follow motor events for mPostStim in frames (used for eventType 3)
bopts.framesPerTrial = vid_end; % nr. of frames per trial
bopts.folds = 10; %nr of folds for cross-validation


% bmat = table2array(events(:,~strcmp(events.Properties.VariableNames, 'Video End')));


% use only first frame times
bmat = zeros(size(events,1), sum(~strcmp(events.Properties.VariableNames, 'Video End')));
count = 1;
for i = 1:size(b_idx,2)
    if strcmp(events.Properties.VariableNames{i}, 'Video End')
        continue
    elseif contains(events.Properties.VariableNames{i}, 'Drop')
        continue
        bmat(:,i) = table2array(events(:,i));        
    else
        behaviors(:, count) = table2array(events(:,i));
        bmat(b_idx{i}(:,1),i) = 1;
        count = count + 1;
    end
end




% consolidate Asymmetric, Bilateral, and Elliptical movements since those
% appear similar in the behavior clustering (louvain) and large bilateral events
% tend to be extremely rare, so this makes it easier for cross-validation
% bilat_vars = {'Elliptical', 'Right Asymmetric', 'Left Asymmetric', 'Bilateral'};
% bilat_idx = [];
% for i = 1 :length(bilat_vars)
%     bilat_idx = [bilat_idx, find(strcmp(events.Properties.VariableNames, bilat_vars{i}))];
% end
% bilat = any(bmat(:, bilat_idx),2);
% bmat(:, bilat_idx) = [];
% bmat = [bmat, bilat];


% get rid of drop
bmat(:,1) = [];

%%
% bmat = table2array(events(:, 3:end-1));

% specify the type of variable for making the design matrix
regLabels = {};
reg_type = [];
for ii = 1:size(events,2)
    if contains(events.Properties.VariableNames{ii}, 'Drop')
        continue
        reg_type = [reg_type 2];
    elseif strcmp(events.Properties.VariableNames{ii}, 'Video End')
        continue
%     elseif any(strcmp(bilat_vars, events.Properties.VariableNames{ii}))
%         continue
    else
        reg_type = [reg_type 3];
    end
    regLabels = [regLabels, events.Properties.VariableNames{ii}];
end

% reg_type = [reg_type, 3]; % add another motor term for the consolidated bilateral movements
% regLabels = [regLabels, {'Bilateral'}];

%% load audio and valve open information

test = readmatrix([event_path, filesep, getAllFiles(event_path, '_trim.txt')]);
audio = test(1:vid_end,5);
valve = test(1:vid_end,6);

flrthresh = flrv>mean(flrv) + std(flrv);
fllthresh = fllv>mean(fllv) + std(fllv);

behaviors = [behaviors, flrthresh, fllthresh];

% get time on for all movement and audio/valve data
audio = get_on_time(audio);
valve = get_on_time(valve);
flrthresh = get_on_time(flrthresh);
flllthresh = get_on_time(fllthresh);

%% 
% concatenate FL movements as binary variables
bmat = [bmat flrthresh fllthresh];
reg_type = [reg_type 3 3];
regLabels = [regLabels, 'FLR', 'FLL'];

%%
% % concatenate 2 STIMULUS variable types for audio tone and valve opening
% bmat = [audio valve bmat(1:vid_end, :)];
% reg_type = [2 2 reg_type];
% regLabels = ['audio', 'valve', regLabels];

%%

% bmat2 = reshape(bmat, [size(bmat,1), 1, size(bmat,2)]);
[dMat, regIdx] = makeDesignMatrix(bmat, reg_type, bopts);
fullR = dMat;

% regLabels = events.Properties.VariableNames(2:end-1);


% [~, activity_timer] = start_timer(any(behaviors,2), 2, fs);
% activity_timer = activity_timer ./ fs;
% fullR = [dMat, activity_timer];
% regIdx = [regIdx; max(regIdx)+1];
% regLabels = [regLabels, 'Timer'];


% add forelimb movements to the design matrix as continuous variables
% fullR = [dMat flrv fllv];
% regIdx = [regIdx; max(regIdx)+1; max(regIdx)+2]; %regressor index



% regLabels = [regLabels, 'FLR', 'FLL'];
% regLabels = [regLabels, 'rFLv', 'lFLv'];


%% run QR and check for rank-defficiency. This will show whether a given regressor is highly collinear with other regressors in the design matrix.
% The resulting plot ranges from 0 to 1 for each regressor, with 1 being
% fully orthogonal to all preceeding regressors in the matrix and 0 being
% fully redundant. Having fully redundant regressors in the matrix will
% break the model, so in this example those regressors are removed. In
% practice, you should understand where the redundancy is coming from and
% change your model design to avoid it in the first place!

rejIdx = false(1,size(fullR,2));
[~, fullQRR] = qr(bsxfun(@rdivide,fullR,sqrt(sum(fullR.^2))),0); %orthogonalize normalized design matrix
figure; plot(abs(diag(fullQRR)),'linewidth',2); ylim([0 1.1]); title('Regressor orthogonality'); drawnow; %this shows how orthogonal individual regressors are to the rest of the matrix
axis square; ylabel('Norm. vector angle'); xlabel('Regressors');
if sum(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1))) < size(fullR,2) %check if design matrix is full rank
    temp = ~(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1)));
    fprintf('Design matrix is rank-defficient. Removing %d/%d additional regressors.\n', sum(temp), sum(~rejIdx));
    rejIdx(~rejIdx) = temp; %reject regressors that cause rank-defficint matrix
end

%% fit model to imaging data
disp("Running ridgeMML cell")
[ridgeVals, dimBeta] = ridgeMML(Nresample', fullR, true); %get ridge penalties and beta weights.
% save([fPath 'dimBeta.mat'], 'dimBeta', 'ridgeVals'); %save beta kernels

%reconstruct imaging data and compute R^2
Vm = (fullR * dimBeta)';

%%
figure


pos_mod = cell(length(regLabels));
ned_mod = cell(length(regLabels));

for j =1:length(regLabels)
    cVar = regLabels(j);%{'Right Asymmetric'};
    clear cIdx
    for i = 1:length(cVar)
        cIdx(i,:) = regIdx == find(ismember(regLabels,cVar{i}));
    end
    cIdx = any(cIdx, 1);
    
    % test = fullR(:,cIdx);
    test = dimBeta(cIdx,:)';
    [~, I] = sort(mean(test(:,0.5*fs:end),2));
    subplot(4, 4, j)
    imagesc(test(I,:))
    colorbar
    caxis([-1 1])
    colormap(bluewhitered())
    try if reg_type(j) == 3 
            vline(0.5*fs, 'k-')
            
            baseline = test(:,1:0.5*fs);
            mu0 = mean(baseline,2);
            std0 = std(baseline,[],2);

%             mux = mean(test(I,0.5*fs+1:end),2);
            response_snr = max(test(:,0.5*fs+1:end), [], 2) ./ mu0;

        end
    catch
    end
    title(regLabels{j})
end

%%

pos_mod = response_snr>mean(response_snr) + 2*std(response_snr);

%%

figure, plot(test(I(pos_mod),:))



%%



figure
for j = 4%1:length(regLabels)
%     cVar = {'Right'};
    cVar = regLabels(j);
    clear cIdx
    for i = 1:length(cVar)
        cIdx(i,:) = regIdx == find(ismember(regLabels,cVar{i}));
    end
    cIdx = any(cIdx, 1);
    
    % test = fullR(:,cIdx);
    test = dimBeta(cIdx,:)';
    [~, I] = sort(mean(test(:,0.5*fs:end),2));
    subplot(4, 4, j)
    imagesc(test(I,:))
    colorbar
%     caxis([-2 4])
    colormap(bluewhitered())
    try if reg_type(j) == 3, vline(0.5*fs, 'k-'), end
    catch
    end
    title(regLabels{j})
end


%%

figure, 
subplot(1,2,1), 
plot(xt(test,fs,2)-0.5, test(I(end-100:end),:)', 'color', [0 0 0 0.25])
hold on
plot(xt(test,fs,2)-0.5, mean(test(I(end-100:end),:)), 'color', 'k', 'LineWidth', 2)
axis([-0.5 2.5 -1 5])
vline(0, 'r-')

subplot(1,2,2), 
plot(xt(test,fs,2)-0.5, test(I(1:100),:)', 'color', [0 0 0 0.25]),
hold on
plot(xt(test,fs,2)-0.5, mean(test(I(1:100),:)), 'color', 'k', 'LineWidth', 2)
% axis tight
axis([-0.5 2 -1 5])
vline(0, 'r-')
%%


%reconstruct imaging data and compute R^2
vv = (fullR(:, cIdx) * dimBeta(cIdx, :))';

%%

figure, hold on
for i = 1:20
    plot(xt(Nresample, fs, 2), Nresample(i,:)-i*5, 'k', 'LineWidth', 1) 
    plot(xt(Nresample, fs, 2), Vm(i,:)-i*5)
end


%%

figure, imagesc(dimBeta)
caxis([-0.05 0.05]), colormap(bluewhitered())
hold on
hline(find(diff(regIdx))+0.5, 'w-')

%%
disp('Running reduced models')
reducedMat = [];
for i = 1:length(regLabels)
    disp(['Starting iteration ', num2str(i)])
    reduced = fullR;
    cIdx = regIdx == i;
    if ~any(cIdx)
        disp(['No ', regLabels{i}, '  events found. Skipping...'])
        continue
    end
    reduced(:, cIdx) = reduced(randperm(size(reduced, 1)), cIdx);

%     [ridgeValsR, dimBetaR] = ridgeMML(Vbrain', reduced, true); %get ridge penalties and beta weights.

    [Vreduced{i}, reducedBeta{i}, reducedR, reducedIdx, reducedRidge, reducedLabels] = crossValModel(reduced, Nresample, regLabels, regIdx, regLabels, bopts.folds);
%     reducedMat(:, i) = modelCorr(Vmaster, Vreduced{i}, Umaster) .^2; %compute explained variance
end

%%

%%

% [~,I] = sort(mysse);

figure, hold on
for i = 1:20
    plot(xt(Nresample, fs, 2), Nresample(I(i),:)-i*4, 'k', 'LineWidth', 1) 
    plot(xt(Nresample, fs, 2), Vm(I(i),:)-i*4)
end


%%
corrMat = modelCorr(Vc,Vm,U) .^2; %compute explained variance
corrMat = arrayShrink(corrMat,mask,'split'); %recreate full frame
corrMat = alignAllenTransIm(corrMat,opts.transParams); %align to allen atlas
corrMat = corrMat(:, 1:size(allenMask,2));

%% run cross-validation
%full model - this will take a moment
disp('Running ridge regression with 10-fold cross-validation')
[Vfull, fullBeta, ~, fullIdx, fullRidge, fullLabels] = crossValModel(fullR, Vbrain, regLabels, regIdx, regLabels, bopts.folds);
save([fPath, char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_cvFull.mat'], 'Vfull', 'fullBeta', 'fullR', 'fullIdx', 'fullRidge', 'fullLabels', '-v7.3'); %save some results

fullMat = modelCorr(Vbrain,Vfull,Ubrain) .^2; %compute explained variance
    
    %% run reduced models for unique contribution
    if size(event_type,1) <= 1
        fig_flag = 0;
        disp('Not enough unique event types to run reduce models. Skipping...')
    else
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
    end

%%

cVar = {'Right'};

for j = 1:length(regLabels)
    cVar = regLabels(j);
    clear cIdx
    for i = 1:length(cVar)
        cIdx(i,:) = regIdx == find(ismember(regLabels,cVar{i}));
    end
    cIdx = any(cIdx, 1);
    
    % test = fullR(:,cIdx);
%     test = mean(catcell(3, reducedBeta{j}),3)';
    test = reducedBeta{j}{1}';
    test = test(:, cIdx);
%     test = dimBeta(cIdx,:)';
    [~, I] = sort(mean(test(:,0.5*fs:end),2));
    subplot(4, 4, j)
    imagesc(test(I,:))
    colorbar
%     caxis([-2 4])
    colormap(bluewhitered())
    try if reg_type(j) == 3, vline(0.5*fs, 'k-'), end
    catch
    end
    title(regLabels{j})
end
%%

meanBetas = mean(catcell(3, reducedBeta{2}),3)';
% [~, I] = sort(mean(meanBetas(:, )))

figure, imagesc(meanBetas)





function new_dat = get_on_time(data)
new_dat = zeros(size(data));
d = diff(data);
t_on = find(d>0)+ 1;

new_dat(t_on) = 1;

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