clc, clear
fileID = fopen('expt3_datalist.txt','r');
addpath('C:\Users\user\Documents\Nick\grooming\utils')

addpath('C:\Users\user\Documents\Nick\ridgeModel');
addpath('C:\Users\user\Documents\Nick\ridgeModel\widefield')
addpath('C:\Users\user\Documents\Nick\ridgeModel\smallStuff') 
addpath('C:\Users\user\Documents\Nick\grooming\utils')

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
data_list = data_list{1};

ger2_idx = 1:11;
hyl3_idx = 52:56;
ecr2_idx = 57:59;


%%
clc

expts_to_analyze = [ger2_idx(1:6)];
for i = 11:length(expts_to_analyze)
    fpath = data_list{expts_to_analyze(i)};
    [events, b_idx, b_tab] = read_boris([fpath, filesep, getAllFiles(fpath, 'events.tsv')]);
    load([fpath, filesep, getAllFiles(fpath, 'cam0_svd.mat')])
    dlc_speed = readmatrix([fpath, filesep, getAllFiles(fpath, '_speed.csv')]);
    ini = ini2struct([fpath, filesep, getAllFiles(fpath, '.ini')]);
    try ini = ini.sentech_give_rewards;
    catch 
        ini = ini.sentech_dlc_live;
    end
    fs = str2double(ini.framerate);

    % consolidate limb speeds from all angles
    fl_l = sqrt(dlc_speed(:,1).^2 + dlc_speed(:,2).^2);
    fl_r = sqrt(dlc_speed(:,3).^2 + dlc_speed(:,4).^2);
    hl_l = dlc_speed(:,5);
    hl_r = dlc_speed(:,6);
    snout = sqrt(dlc_speed(:,7).^2 + dlc_speed(:,8).^2 + dlc_speed(:,9).^2);
    tailbase = sqrt(dlc_speed(:,10).^2 + dlc_speed(:,11).^2);
    
    all_movement = sqrt(fl_l.^2 + fl_r.^2 + hl_l.^2 + hl_r.^2 + snout.^2 + tailbase.^2);
    mvt = resample(all_movement, size(V,2), length(all_movement));


    % filter brain data
    [b, a] = butter(1, [0.01 10]/(fs/2));
    Vbrain = s*V;
    Vbrain = filtfilt(b,a, Vbrain')';


    disp('Building design matrix')
    bopts.frameRate = fs;
    bopts.sPostTime=round(fs*1);
    bopts.mPreTime = ceil(0.5 * fs);  % precede motor events to capture preparatory activity in frames (used for eventType 3)
    bopts.mPostTime = ceil(2 * fs);   % follow motor events for mPostStim in frames (used for eventType 3)
    bopts.framesPerTrial = size(Vbrain,2); % nr. of frames per trial
    bopts.folds = 1;
    
    
    
    mvtOn = get_on_time(mvt>(mean(mvt)+std(mvt)));
    lick_start = zeros(size(mvt));
    lick_start(b_idx{1}(:,1)) = 1;
end


%% consolidate limb speeds from all angles

fl_l = sqrt(dlc_speed(:,1).^2 + dlc_speed(:,2).^2);
fl_r = sqrt(dlc_speed(:,3).^2 + dlc_speed(:,4).^2);
hl_l = dlc_speed(:,5);
hl_r = dlc_speed(:,6);
snout = sqrt(dlc_speed(:,7).^2 + dlc_speed(:,8).^2 + dlc_speed(:,9).^2);
tailbase = sqrt(dlc_speed(:,10).^2 + dlc_speed(:,11).^2);

all_movement = sqrt(fl_l.^2 + fl_r.^2 + hl_l.^2 + hl_r.^2 + snout.^2 + tailbase.^2);
mvt = resample(all_movement, size(V,2), length(all_movement));

%%

ll = zeros(1,size(V,2));
for i = 1:size(b_idx{1},1)
ll(b_idx{1}(i,1):b_idx{1}(i,2)) = 1;
end

%%

[b, a] = butter(2, 0.01/(fs/2), 'high');


pc = 1;
negate = 1;
test = filtfilt(b, a, V(pc,1:end-1));
% test = V(2,1:end-1);

t = xt(test, fs, 2);

pcmap = imrotate(reshape(U(:,pc), [128 128]),-90);

if negate
    test = -test;
    pcmap = -pcmap;
end


figure, 
subplot(2,5,[1,2,6,7])
imagesc(pcmap)
colorbar

subplot(2,5,3:5)
plot(t, test), hold on, patchplot(t(b_idx{1}), [min(test) max(test)], 'm')
pp = 20e-4*(mvt>mean(mvt));
pp(pp==0) = nan;
plot(t, pp(1:end-1), 'k', 'LineWidth', 2)
% plot(t)
axis tight
xlabel('Time (s)')
ylabel('Filtered fluorescence')


[~,f,tt, ps] = spectrogram(test,150,50,150,fs,'yaxis');
subplot(2,5, 8:10)
psdb = 10*log10(ps);
pcolor(tt,f,psdb), shading interp
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cc = colorbar;
cc.Label.String = 'Power (dB)';
caxis([prctile(psdb,50,'all'), max(psdb(:))])



%%


figure, plot(zscore(test)-2), hold on, plot(zscore(all_movement))

%% Hypothesis:
%  Continuous licking behavior inhibits dorsal cortex
%  Method:
%       Segment trial data into drinking and non-drinking periods
%       Run regression on brain and behavior data for each segment
%           Keep U the same between both
%           Compare Betas magnitude and shape


% Aggregation window
W = 1;

% Resample behavior so it's the same length as the brain
mvt = resample(all_movement, size(V,2), length(all_movement));


% filter brain data
[b, a] = butter(1, [0.01 10]/(fs/2));
Vbrain = s*V;
Vbrain = filtfilt(b,a, Vbrain')';

% segment brain data based on behavior - create binary lick data
ll = zeros(1,size(V,2));
ll(b_idx{1}) = 1;


% Vbrain_licking = Vbrain(:,logical(ll));
% mvt_licking = mvt(logical(ll));
% 
% Vbrain_notlicking = Vbrain(:, ~logical(ll));
% mvt_notlicking = mvt(~logical(ll));


%%
disp('Building design matrix')
bopts.frameRate = fs;
bopts.sPostTime=round(fs*1);
bopts.mPreTime = ceil(0.5 * fs);  % precede motor events to capture preparatory activity in frames (used for eventType 3)
bopts.mPostTime = ceil(2 * fs);   % follow motor events for mPostStim in frames (used for eventType 3)
bopts.framesPerTrial = size(Vbrain,2); % nr. of frames per trial
bopts.folds = 1;



mvtOn = get_on_time(mvt>(mean(mvt)+std(mvt)));
lick_start = zeros(size(mvt));
lick_start(b_idx{1}(:,1)) = 1;
% cont_mvt = ll(1:length(mvt));

[movMat, regIdx] = makeDesignMatrix([mvtOn, lick_start], [3 3], bopts);
fullR = [movMat];

% fullR = mvt_licking;
% regIdx = 1;

regLabels = {'Movement', 'Lick'};
disp('Running ridge regression with 10-fold cross-validation')
[Vfull, fullBeta, ~, fullIdx, fullRidge, fullLabels] = crossValModel(fullR, Vbrain, regLabels, regIdx, regLabels, bopts.folds);
% save([fPath, char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_cvFull.mat'], 'Vfull', 'fullBeta', 'fullR', 'fullIdx', 'fullRidge', 'fullLabels', '-v7.3'); %save some results

% fullMat = modelCorr(Vbrain_licking,Vfull,U) .^2; %compute explained variance

%%

% [L, betas, convergenceFailures] = ridgeMML(Vbrain_licking(1:10,:)', mvt_licking)




fullMat = modelCorr(Vbrain,Vfull,U) .^2;


%%

visual = true;
cBetaRight = check_beta('Movement', fullLabels, fullIdx, U, fullBeta{1}, Vfull, [], visual);
% right = movmean(cBetaRight, 6, 3);




function new_dat = get_on_time(data)
new_dat = zeros(size(data));
d = diff(data);
t_on = find(d>0)+ 1;

new_dat(t_on) = 1;

end
