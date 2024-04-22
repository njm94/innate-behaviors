%% Compile behavior data
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

%% test section
f = figure;
f.Position = [100 100 1100 500];
h = 30;

for ii = 1:length(beh)
    for i = 1:size(beh{ii},4)
        if i-h < 1
            a = 1;
        else
            a = i-h;
        end
        subplot(1,2,1)
        imagesc(beh{ii}(:,:,:,i)), 
        hold on,
        plot(paw{ii}(a:i,1), paw{ii}(a:i,2), 'r', 'LineWidth', 2),  axis off
        plot(paw{ii}(a:i,3), paw{ii}(a:i,4), 'b', 'LineWidth', 2)
        
        [I,J] = find(mean(led_behavior_roi,3));
        line([min(J) max(J)], [min(I) min(I)], 'color', [1 1 1], 'LineWidth', 2)
        line([min(J) max(J)], [max(I) max(I)], 'color', [1 1 1], 'LineWidth', 2)
        line([min(J) min(J)], [min(I) max(I)], 'color', [1 1 1], 'LineWidth', 2)
        line([max(J) max(J)], [min(I) max(I)], 'color', [1 1 1], 'LineWidth', 2)
        
        [I,J] = find(mean(whisker_behavior_roi,3));
        line([min(J) max(J)], [min(I) min(I)], 'color', [1 1 1], 'LineWidth', 2)
        line([min(J) max(J)], [max(I) max(I)], 'color', [1 1 1], 'LineWidth', 2)
        line([min(J) min(J)], [min(I) max(I)], 'color', [1 1 1], 'LineWidth', 2)
        line([max(J) max(J)], [min(I) max(I)], 'color', [1 1 1], 'LineWidth', 2)
        hold off
        
        
        text(5, 10, ['Trial: ', num2str(ii)], 'FontSize', 16, 'Color', [1 1 1])
        if whisker_time{ii}(i)
            title('WHISKER')
        elseif led_behavior_time{ii}(i)
            title('LED')
        elseif rewardTime(ii) ~= 1 && abs(i/fs - rewardTime(ii)) < 0.15
            title('REWARD')
        else
            title('')
        end
%         suptitle(num2str(i))

        subplot(1,2,2),
        imagesc(dFF{ii}(:,:,i)), colormap jet, colorbar, caxis([-0.1 0.1]), axis off
        pause(0.01)
    end
end
%%
close all
fs = 60;
f = figure;
f.Position = [100 100 1100 500];
h = 30;
count = 1;
for ii = 1:length(beh)
    for i = 1:size(beh{ii},4)
        if i-h < 1
            a = 1;
        else
            a = i-h;
        end
        subplot(1,2,1)
        imagesc(beh{ii}(:,:,:,i)), 
        hold on,
        plot(paw{ii}(a:i,1), paw{ii}(a:i,2), 'r', 'LineWidth', 2),  axis off
        plot(paw{ii}(a:i,3), paw{ii}(a:i,4), 'b', 'LineWidth', 2)
        
        [I,J] = find(mean(led_behavior_roi,3));
        line([min(J) max(J)], [min(I) min(I)], 'color', [1 1 1], 'LineWidth', 2)
        line([min(J) max(J)], [max(I) max(I)], 'color', [1 1 1], 'LineWidth', 2)
        line([min(J) min(J)], [min(I) max(I)], 'color', [1 1 1], 'LineWidth', 2)
        line([max(J) max(J)], [min(I) max(I)], 'color', [1 1 1], 'LineWidth', 2)
        
        [I,J] = find(mean(whisker_behavior_roi,3));
        line([min(J) max(J)], [min(I) min(I)], 'color', [1 1 1], 'LineWidth', 2)
        line([min(J) max(J)], [max(I) max(I)], 'color', [1 1 1], 'LineWidth', 2)
        line([min(J) min(J)], [min(I) max(I)], 'color', [1 1 1], 'LineWidth', 2)
        line([max(J) max(J)], [min(I) max(I)], 'color', [1 1 1], 'LineWidth', 2)
        hold off
        
        
        text(5, 10, ['Trial: ', num2str(ii)], 'FontSize', 16, 'Color', [1 1 1])
        if whisker_time{ii}(i)
            title('WHISKER')
        elseif led_behavior_time{ii}(i)
            title('LED')
        elseif rewardTime(ii) ~= 1 && abs(i/fs - rewardTime(ii)) < 0.15
            title('REWARD')
        else
            title('')
        end
%         suptitle(num2str(i))

        subplot(1,2,2),
        imagesc(dFF{ii}(:,:,i)), colormap jet, caxis([-0.1 0.1]), axis off
        c=colorbar; c.Label.String = '\DeltaF/F_0'; 
%         pause(0.01)
        F(count) = getframe(f);
        count = count + 1;
    end
end

%%
% savePath = 'Y:\pankaj\water_reaching_project\data\matlab_outputs';
writerObj = VideoWriter([savePath, filesep, ...
            mouseID{1}, '_', uDays{1}, 'all_trials', '.avi']);
writerObj.FrameRate = 60;

  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
disp('Saving...')
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
disp('Done')
%% create a structure of length mouseID which contains all data organized
fs = 60;
dur = 10;

for i = 4%:numel(mouseID) % loop over mice   
    
    
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
    for j = 1%:length(uDays) % loop over days

%         if isfile([savePath, filesep, mouseID{i}, '_', uDays{j}, 'data_compile', '.mat'])
%             disp([savePath, filesep, mouseID{i}, '_', uDays{j}, 'data_compile', '.mat already exists. Skipping...'])
%             continue
%         end
        % BBMAP FILE INFO
        [behPath, brPath, trialResult, brFs, behFs, rewardTime, lickTime] = getTrialInfo(mouseData{i}, uDays{j});
        
        rois_defined = false;
        led_stim_is_good = false;

        for k = 1:length(behPath) % loop over trials           
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
            [fluo{k}, ref{k}] = splitstrobe(I);
            fluo_frame = fluo{k}(:,:,10);

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
            [led_behavior_time{k}, led_behavior_roi, led_behavior_signal{k}] = get_stimulus_time(beh{k}, 'LED stim', led_behavior_roi, 1,[230 270], 20);  
            rois_defined = true;
            
%           % find flash artifact in upper left corner of image
%             disp('Finding LED flash artifact');
% %             [led_brain_time{k}, led_brain_roi, led_brain_signal{k}] = get_stimulus_time(fluo{k}, 'LED flash', led_brain_roi, 1,[230 270], 20);
%             while ~led_stim_is_good            
%                 [led_brain_time{k}, led_brain_roi, led_brain_signal{k}] = get_stimulus_time(fluo{k}, 'LED flash', led_brain_roi, 1,[230 270], 20);
%                 h=figure;
%                 plot(led_brain_signal{k})
%                 answer = questdlg('is it good?');
%                 switch answer
%                     case 'Yes'
%                         led_stim_is_good = true;
%                     otherwise                        
%                         led_brain_roi = [];
%                         led_stim_is_good = false;
%                 end
%                 close(h)
%             end
            
            % calculate dFF based on baseline (before water drop) assuming
            % framerate is 60fps
            if rewardTime(k) ~= -1
                baseline_end = round(rewardTime(k)*fs);
            else
                baseline_end = size(fluo{k},3)-5;
            end
            FdFF = (single(fluo{k}) - mean(fluo{k}(:,:,1:baseline_end),3))./mean(fluo{k}(:,:,1:baseline_end),3);
            RdFF = (single(ref{k}) - mean(ref{k}(:,:,1:baseline_end),3))./mean(ref{k}(:,:,1:baseline_end),3);    
            dFF{k} = (FdFF - RdFF);
            
           
            % remove dark frame artifact at end of trial. 5 frames should
            % be enough
            dFF{k} = dFF{k}(:,:,1:end-5);
            min_size = min([size(dFF{k},3) size(beh{k},4)]);
            dFF{k} = dFF{k}(:,:,1:min_size);
            beh{k} = beh{k}(:,:,:,1:min_size);
            led_behavior_time{k} = led_behavior_time{k}(1:min_size);
%             led_brain_time{k} = led_brain_time{k}(1:min_size);
            whisker_time{k} = whisker_time{k}(1:min_size);
            paw{k} = paw{k}(1:min_size,:);
         
        end

        disp('Saving ...')
        tic
%         save([savePath, filesep, ...
%             mouseID{i}, '_', uDays{j}, 'data_compile', '.mat'], ...
%             'group', 'behPath', 'brPath', 'trialResult', 'brFs', ...
%             'behFs', 'rewardTime', 'lickTime', 'beh', 'paw', ...
%             'whisker_time', 'led_behavior_time', 'dFF', ...
%             'led_brain_time', 'fluo', 'ref', 'fluo_frame', ...
%             'led_brain_signal', 'led_behavior_signal', ...
%             'led_brain_roi', 'led_behavior_roi', 'whisker_behavior_roi', ...
%             '-v7.3')
        save([savePath, filesep, ...
            mouseID{i}, '_', uDays{j}, 'data_compile', '.mat'], ...
            'group', 'behPath', 'brPath', 'trialResult', 'brFs', ...
            'behFs', 'rewardTime', 'lickTime', 'beh', 'paw', ...
            'whisker_time', 'led_behavior_time', 'dFF', ...
            'fluo', 'ref', 'fluo_frame', 'led_behavior_signal', ...
            'led_behavior_roi', 'whisker_behavior_roi', '-v7.3')
        toc
    end
    
        
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




% return structure of length uniqueDays containing relevant trial info
% trials(length(uDays)) = struct();
% for i = 1:length(uDays)
%     idx = contains(behPath, [uDays{i},'_']);
%     trials(i).behPath = behPath(idx);
%     trials(i).brPath = brPath(idx);
%     trials(i).trialResult = trialResult(idx);
%     trials(i).brFs = cell2mat(brFs(idx));
%     trials(i).behFs = cell2mat(behFs(idx));
%     trials(i).rewardTime = cell2mat(rewardTime(idx));
%     trials(i).brFrames = cell2mat(brFrames(idx));
%     trials(i).lickTime = cell2mat(lickTime(idx));
% end

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
    bfiles = bfiles(contains(bfiles, '.csv') & contains(bfiles, bfileID{:}));
    tab = readtable([bpath, bfiles{:}]);

    l_paw_x = table2array(tab(1:end, 5));
    l_paw_y = table2array(tab(1:end, 6));

    r_paw_x = table2array(tab(1:end, 47));
    r_paw_y = table2array(tab(1:end, 48));
    
%     dlx = diff(l_paw_x);
%     dly = diff(l_paw_y);
%     lv = sqrt(dlx.^2 + dly.^2);
%     drx = diff(r_paw_x);
%     dry = diff(r_paw_y);
%     rv = sqrt(drx.^2 + dry.^2);
%     
%     paw.left.pos = [l_paw_x, l_paw_y];
%     paw.left.vel = lv;
%     paw.right.pos = [r_paw_x, r_paw_y];
%     paw.right.vel = rv;
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
function [artifact, stimCenter] = get_artifacts(fluo, ref)
    tmp = movstd(diff(fluo),3);
    thresh = mean(tmp) + 5*std(tmp);
    artifact1 = tmp>thresh;

    tmp = movstd(diff(ref),3);
    thresh = mean(tmp) + 5*std(tmp);
    artifact2 = tmp>thresh;

    artifact = artifact1 | artifact2;
    artifact = logical(cat(1, artifact, 0));
    [~, stimCenter] = findpeaks(smoothdata(artifact,'gaussian',10));
end
function I = load_n_frames(tiff_file, n)

info = imfinfo(tiff_file);
tampon = imread(tiff_file,'Index',1);
F = length(info);

if nargin < 2 || isempty(n)
    n = F;
elseif n > F
    error('Trying to load more frames than exists in the file')
end


I = zeros(size(tampon,1),size(tampon,2),n,'uint16');
I(:,:,1) = tampon(:,:,1);
tic
wait_bar = waitbar(0,['Loading ',tiff_file]);
ind = 0;
for i = 2:n
    if ind == 0, waitbar(i/n, wait_bar); end
    ind = ind + 1; if ind == 100, ind = 0; end
    tampon = imread(tiff_file,'Index',i,'Info',info);
    I(:,:,i) = tampon(:,:,1);
end
close(wait_bar);
temps = num2str(round(10*toc)/10);
disp([tiff_file ' open in ' num2str(temps) 's'])
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