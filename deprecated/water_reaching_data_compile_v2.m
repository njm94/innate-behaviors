%% Compile behavior data
%
% Modification on June 26 2024. Do PCA on brain data and filter and zscore
% in the same way to match grooming data. Save the fluorescence alone as
% well as the hemodynamic corrected data.
% 
% This code separates fluorescence and reflectance images, reads the brain
% behavior map file, and reads the deeplabcut track file. The code does not
% do LED artifact correction, removal of dark frames at the end of the
% trial, calculate dFF, or any quality control.
% 
% FROM EXPERIMENT INFORMATION
%   group - group number (one of [1, 2, 3, 4] )
%
% FROM BRAIN BEHAVIOR MAP FILE
%   trial_result
%       (cell array of length num_trials)
%   reward_time - time when drop appears, expressed in seconds
%       (array of length num_trials)
%   lick_time - time when touch is detected, expressed in seconds
%       (array of length num_trials)
%       (sometimes touch detector does not work!! verify with trial_result)
%   brain_est_fps - estimated frame rate for brain imaging
%       (These values are based on pi and may not be trustworthy)
%   beh_est_fps - estimated frame rate for behavior imaging
%       (These values are based on pi and may not be trustworthy)
%   brain_path - path to original brain data
%   beh_path - path to original behavior data
%
% FROM BEHAVIOR VIDEO DATA
%   beh - behavior video
%       (cell array of length num_trials, each element is size
%           320 x 320 x 3 x num_behavior_frames)
%   whisker_time - whisker stimulus time, expressed in seconds
%       (cell array of length num_trials, each element is logical array of
%       size [num_behavior_frames 1], empty if no whisker stimulus)
%   LED_behavior_time - LED stimulus time, expressed in frames
%       (array of length num_trials)
%
% FROM DEEPLABCUT OUTPUTS
%   paw - paw kinematics
%       (cell array of length num_trials, each element is size
%           num_behavior_frames x 4, where each column is:
%           L_paw_x, L_paw_y, R_paw_x, R_paw_y)
%
% FROM BRAIN IMAGING
%   fluo - fluorescence data for brain imaging
%       (cell array of length num_trials, each element is size 
%           128 x 128 x num_brain_frames)
%   ref - reflectance data for brain imaging 
%       (cell array of length num_trials, each element is size 
%           128 x 128 x num_brain_frames)
%   fluo_frame - single frame from fluorescence data for each trial
%       (array of size 128 x 128 x num_trials)



% EXPERIMENT INFORMATION
% GROUP 1
%   HR1   Stroke
%   HR2   Stroke
%   IK1  Stroke
%   IK2  Stroke
% GROUP 2
%   B1   Stroke No brain
%   B2   Stroke No brain
%   MB2   Stroke No brain
% GROUP 3
%   IJ1  Sham   No brain
%   IJ2  Stroke
%   IJ3  Sham
%   IO1  Stroke
%   IQ1  Stroke
%   IQ2  Sham
%   IQ3  Sham   No brain

clc, clear
addpath('Y:\pankaj\water_reaching_project\widefield_imaging')

mapFile = 'Y:\pankaj\\water_reaching_project\data\brain_behavior_map_all_20221126.xlsx';
savePath = 'Y:\pankaj\water_reaching_project\data\matlab_outputs';
[~, mouseID] = xlsfinfo(mapFile);

% Save trial information into a cell array where each cell corresponds to a
% single mouse
mouseData = cell(numel(mouseID), 1);
for i = 1:numel(mouseID)
    mouseData{i} = readtable(mapFile, 'Sheet', mouseID{i});
end

%% create a structure of length mouseID which contains all data organized
fs = 60;
dur = 10;
num_PCs = 500;


[b, a] = butter(2, 0.01/(fs/2), 'high');

for i = 4%:4%:13%:numel(mouseID) % loop over mice       
    % FROM EXPERIMENT INFO. IJ1 and IQ3 are placed in group2 because of
    % no-brain data
    switch mouseID{i}
        case {'HR1', 'HR2', 'IK1', 'IK2'}
            group = 1;
        case {'B1', 'B2', 'MB2', 'IJ1', 'IQ3'}
            group = 2;
        case {'IJ2', 'IJ3', 'IO1', 'IQ1', 'IQ2'}
            group = 3;
    end

    % Group 2 did not contain any brain data. Skip
    if group == 2, continue; end    
    
    uDays = getUniqueDays(mouseData{i});
    whisker_time = cell(length(uDays),1);
    paw = cell(length(uDays),1);
    led_behavior_time = cell(length(uDays),1);
    led_behavior_signal = cell(length(uDays),1);
    beh = cell(length(uDays),1);
    dFF = cell(length(uDays),1);
    led_brain_time = cell(length(uDays),1);
    led_brain_signal = cell(length(uDays),1);
    fluo = cell(length(uDays),1);
    ref = cell(length(uDays),1);


    % BBMAP FILE INFO
    [behPath, brPath, trialResult, brFs, behFs, rewardTime, lickTime] = getTrialInfo(mouseData{i}, uDays{1});
    
    rois_defined = false;
    led_stim_is_good = false;

    % We are only interested in day 1, as the rest of the days involve
    % stroke. Loop over trials in day 1
    for k = 1:length(behPath)     
        if ~strcmpi(trialResult{k}, 'success'), continue; end
        % FROM DEEPLABCUT OUTPUT
        if contains(behPath{k}, '/media/pankaj/teamshare/')
            behPath{k} = strrep(behPath{k}, '/media/pankaj/teamshare', 'Y:');
            behPath{k} = strrep(behPath{k}, '/', '\');
        end
        
        path_separators = strfind(behPath{k}, '\');
        behPath{k} = [behPath{k}(1:path_separators(2)), 'to_merge', behPath{k}(path_separators(2):end)];
        paw{k} = read_DLC_file(behPath{k});
        
        % behavior video
        v = VideoReader([behPath{k}(1:end-4) 'avi']);
        beh{k} = read(v, [1 Inf]);                      
        
        if contains(brPath{k}, '/media/deeplabcut/data2/')
            brPath{k} = strrep(brPath{k}, '/media/deeplabcut/data2', 'Y:\pankaj');
            brPath{k} = strrep(brPath{k}, '/', '\');
        end

        if contains(brPath{k}, '/media/pankaj/teamshare')
            brPath{k} = strrep(brPath{k}, '/media/pankaj/teamshare', 'Y:');
            brPath{k} = strrep(brPath{k}, '/', '\');
        end
        
        I = loadtiff(brPath{k});
        [ff, rr] = splitstrobe(I);

        % remove few frames from the end because sometimes there is a light
        % off artifact
        ff = ff(:,:,1:end-3);
        rr = rr(:,:,1:end-3);


        fluo_frame = ff(:,:,10);

        % do SVD on brain data
        [Uf, sf, Vf] = movie_svd(ff, num_PCs);
        [Ur, sr, Vr] = movie_svd(rr, num_PCs);
        Vf = filtfilt(b, a, (sf*Vf')')';
        Vr = filtfilt(b, a, (sr*Vr')')';

        fluo{k} = reshape(Uf * Vf, [size(ff,1), size(ff,2), size(ff,3)]);

        ref{k} = reshape(Ur * Vr, [size(rr,1), size(rr,2), size(rr,3)]);

        % stimulus timestamps
        if ~rois_defined
            whisker_behavior_roi = [];
            led_behavior_roi = [];
            led_brain_roi=[];
        end
%             
        if rewardTime(k) ~= -1 && group == 1
            [whisker_time{k}, whisker_behavior_roi] = get_stimulus_time(beh{k}, 'whisker stim', whisker_behavior_roi, 3.5, [100 200], 40);
        else
            whisker_time{k} = false(size(beh{k},4),1);
        end
        rois_defined = true;     
        
        % calculate dFF based on baseline (before water drop) assuming
        % framerate is 60fps
        if rewardTime(k) ~= -1
            baseline_end = round(rewardTime(k)*fs);
        else
            baseline_end = size(fluo{k},3);
        end

        % since data is filtered, subtract min to bring mean away from 0 
        % before calculating dFF
        FdFF = fluo{k} - min(fluo{k}(:));
        FdFF = (FdFF ./ mean(FdFF(:,:,1:baseline_end), 3)) - 1;
        
        RdFF = ref{k} - min(ref{k}(:));
        RdFF = (RdFF ./ mean(RdFF(:,:,1:baseline_end), 3)) - 1;

        dFF{k} = zscore((FdFF - RdFF), [], 3);
        FdFF = zscore(FdFF,[],3);

        min_size = min([size(dFF{k},3) size(beh{k},4)]);
        dFF{k} = dFF{k}(:,:,1:min_size);
        beh{k} = beh{k}(:,:,:,1:min_size);
        dff_fluo{k} = FdFF(:,:,1:min_size);
%             led_behavior_time{k} = led_behavior_time{k}(1:min_size);
%             led_brain_time{k} = led_brain_time{k}(1:min_size);
%             whisker_time{k} = whisker_time{k}(1:min_size);
        paw{k} = paw{k}(1:min_size,:);
     
    end

    disp('Saving ...')
    tic

    save([savePath, filesep, ...
        mouseID{i}, '_', uDays{1}, 'data_compile', '.mat'], ...
        'group', 'behPath', 'brPath', 'trialResult', 'brFs', ...
        'behFs', 'rewardTime', 'lickTime', 'beh', 'paw', ...
        'whisker_time', 'led_behavior_time', 'dFF', 'dff_fluo', ...
        'fluo', 'ref', 'fluo_frame', 'led_behavior_signal', ...
        'led_behavior_roi', 'whisker_behavior_roi', '-v7.3')
    toc
    
        
end


%%
function uDays = getUniqueDays(BBMapTable)
% get brain and behavior pats
behPath = table2cell(BBMapTable(:,1));
brPath = table2cell(BBMapTable(:,2));

% pull out the date information
behDays = regexp(behPath, '20\d\d_\d\d?_\d\d?', 'match');
brDays = regexp(brPath, '20\d\d_\d\d?_\d\d?', 'match');

% organize into cell array
behDays = vertcat(behDays{:});
brDays = vertcat(brDays{:});

% get unique values
uBehDays = unique(behDays, 'stable');
uBrDays = unique(brDays, 'stable');


% check that same days exist in brain behavior
assert(length(uBehDays) == length(uBrDays), ...
    'Unequal number of unique days in brain and behavior dataset.');


% make sure days match
for i = 1:length(uBehDays)
    assert(strcmp(uBehDays{i}, uBrDays{i}), ...
        ['Unique day at index ', num2str(i), ...
        ' of behavior data does not agree with brain data.']);        
end


% return cell array of unique days
uDays = uBehDays;

end
function [behPath, brPath, trialResult, brFs, behFs, rewardTime, lickTime] = getTrialInfo(BBMapTable, uDay)
all_behavior_paths = table2cell( BBMapTable(:, 1)  );
idx = contains(all_behavior_paths, [uDay,'_']);

% get paths and trial results
behPath      =   table2cell( BBMapTable(:, 1)  );
brPath       =   table2cell( BBMapTable(:, 2)  );
trialResult  =   table2cell( BBMapTable(:, 3)  );
brFs         =   table2cell( BBMapTable(:, 10) );
behFs        =   table2cell( BBMapTable(:, 6)  );
rewardTime   =   table2cell( BBMapTable(:, 7)  );
lickTime     =   table2cell( BBMapTable(:, 8)  );
% brFrames     =   table2cell( BBMapTable(:, 9)  );

behPath = behPath(idx);
brPath = brPath(idx);
trialResult = trialResult(idx);
brFs = cell2mat(brFs(idx));
behFs = cell2mat(behFs(idx));
rewardTime = cell2mat(rewardTime(idx));
% brFrames = cell2mat(brFrames(idx));
lickTime = cell2mat(lickTime(idx));


end
function [rewardTrialFlag, failedTrialFlag, ignoreFlag] = setFlags(trialResultString)

switch lower(trialResultString)
    case 'ignore'
        ignoreFlag = true;
        rewardTrialFlag = false;
        failedTrialFlag = false;
    case 'donot consider no water'
        ignoreFlag = true;
        rewardTrialFlag = false;
        failedTrialFlag = false;
    case 'no-reward'
        ignoreFlag = false;
        rewardTrialFlag = false;
        failedTrialFlag = false;
    case 'fail'
        ignoreFlag = false;
        rewardTrialFlag = true;
        failedTrialFlag = true;
    case 'miss'
        ignoreFlag = false;
        rewardTrialFlag = true;
        failedTrialFlag = true;
    case 'grooming'
        ignoreFlag = false;
        rewardTrialFlag = true;
        failedTrialFlag = true;
    case 'no move'
        ignoreFlag = false;
        rewardTrialFlag = true;
        failedTrialFlag = true;
    otherwise
        ignoreFlag = false;
        rewardTrialFlag = true;
        failedTrialFlag = false;
end

end    
function paw = read_DLC_file(behPath)
    bfileID = regexp(behPath,'\d{10}','match');
    idx = strfind(behPath, '\');
    bpath = ['Y:', behPath(idx(1):idx(end))];
    bfiles = getAllFiles(bpath);
    bfiles = bfiles(contains(bfiles, '.csv') & contains(bfiles, bfileID{:}) & contains(bfiles, 'DeepCut'));
    tab = readtable([bpath, bfiles{:}]);

    l_paw_x = table2array(tab(1:end, 5));
    l_paw_y = table2array(tab(1:end, 6));

    r_paw_x = table2array(tab(1:end, 47));
    r_paw_y = table2array(tab(1:end, 48));
    
    paw = [l_paw_x, l_paw_y, r_paw_x, r_paw_y];
end
function [stimulus_time, roi, stimulus_signal] = get_stimulus_time(video, stimulus_type, roi, mult, loc, k)
if isempty(roi)
    roi = draw_roi(video, 1, stimulus_type); 
    roi = roi(:,:,:,1);
end
if nargin<4 || isempty(mult), mult=1; end
if nargin<5 || isempty(loc), loc=[1 round(size(video,4)/2)]; end

if length(size(video)) == 4 % behavior video
    stimulus_signal = squeeze(mean(mean(mean(single(video(:,:,:,loc(1):loc(2))).* roi))));
    stim_signal = movmean(abs(squeeze(mean(mean(diff(mean(single(video).* roi,3),1,4))))),k).^3;
elseif length(size(video)) == 3 % brain video
    stimulus_signal = squeeze(mean(mean(single(video(:,:,loc(1):loc(2))).* roi(:,:,1))));
    stim_signal = movmean(abs(squeeze(mean(mean(diff(single(video).* roi(:,:,1),1,3))))),k).^3;
end
% if nargout>2, figure, plot(stimulus_signal); end

thresh = std(stim_signal(1:end-10)); % ignore last few frames in case of brain data where lights turn off

artifact = false(size(stim_signal));
artifact(loc(1):loc(2)) = stim_signal(loc(1):loc(2))>thresh*mult;
idx = find(artifact);
stimulus_time = false(size(stim_signal));
stimulus_time(idx(1):idx(end)) = true;
stimulus_time = [false; stimulus_time];
% figure, plot(stim_signal), hold on, plot(ones(size(stim_signal))*mult*thresh);
end


function varargout = splitstrobe(data,order)
% Split strobed image series into two streams
% Inputs:
%   data  (required - data containing sequentially illuminated frames)
%   order  (optional - frame start (1 or 2) for image series, default = 1)
% Outputs:
%   stream1 (image series consisting of frames 1,3,5,...)
%   stream2 (image series consisting of frames 2,4,6,...)
% Usage: [stream1, stream2] = splitstrobe2(data,order)

if nargin < 2, order = 1; end

N = length(data);
index = order:2:N;
varargout{1} = data(:,:,index);
if nargout > 1
    if order == 2
        opp = 1;
    else
        opp = 2;
    end
    index = opp:2:N;
    varargout{2} = data(:,:,index);
    
    if size(varargout{2},3) > size(varargout{1},3)
        varargout{2} = varargout{2}(:,:,1:end-1);
    elseif size(varargout{1},3) > size(varargout{2},3)
        varargout{1} = varargout{1}(:,:,1:end-1);
    end
end
end