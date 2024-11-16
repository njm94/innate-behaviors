%% behavior analysis scrap
% 


clear, clc
fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
data_list = textscan(fileID, formatSpec);

fs = 90;

thy1_idx = 1:7;
ai94_idx = 8:13;
camk_idx = 14:25;

% create transition matrix
% states = ["Stop", "Elliptical", "Bilateral Asymmetric", "Bilateral", "Unilateral", "Lick"];
states = ["Stop", "Elliptical", "Asymmetric", "Bilateral", "Unilateral"];


spontaneous = [1,2,8,9,14,15,20,21];
evoked = [3:7,10:13,16:19,22:25];

%%

sessions = {'spontaneous', 'evoked'};

for i = 1:2
    tmat = zeros(numel(states));
    event_raster = [];

    for j = eval(sessions{i})%14:length(data_list{1})
        data_dir = data_list{1}{j};
        [mouse_root_dir, exp_date, ~] = fileparts(data_dir);
        [~, mouse_id, ~] = fileparts(mouse_root_dir);
        [expt_root_dir, ~, ~] = fileparts(mouse_root_dir);
       
        snippets_dir = [data_dir, filesep, 'snippets'];
        [snippets, labels] = parse_snippets(snippets_dir);
        num_events = cellfun(@(data) size(data, 1), snippets);
    
        event_duration{j} = cellfun(@(data) diff(data,1,2), snippets, 'UniformOutput', false);
    
        bmat = zeros(1, max(catcell(2, cellfun(@max, snippets, 'UniformOutput', false))));
        for ii = 1:length(snippets)
            for jj = 1:size(snippets{ii},1)
                bmat(snippets{ii}(jj,1):snippets{ii}(jj,2)) = ii;
            end
        end
    
        % consolidate unilateral events
        bmat2 = bmat;
        bmat2(bmat2==3)=2;
        bmat2(bmat2==4)=3;
        bmat2(bmat==5 | bmat==6) = 4;
        bmat2(bmat==7) = 0;% 5;
        bmat2 = bmat2 + 1;
        
        [episodes, idx] = aggregate(bmat, 3);
        

        for ii = 1:size(idx,1)
            tmp = bmat2(idx(ii,1):idx(ii,2));
            event_idx = [1 tmp] > 1;
            event_start = find(diff(event_idx) == 1);
            for jj = 1:length(event_start)
                if jj == 1
                    tmat(1, tmp(event_start(jj))) = tmat(1, tmp(event_start(jj)))+1;
                else
                    tmat(last_event, tmp(event_start(jj))) = tmat(last_event, tmp(event_start(jj))) + 1;
                end
                last_event = tmp(event_start(jj));
            end
            tmat(last_event, 1) = tmat(last_event, 1) + 1;
        end


    end



    figure
    subplot(1,2,1)
    imagesc(tmat)
    xticks(1:7), xticklabels(states)
    yticks(1:7), yticklabels(states)
    ylabel('From')
    xlabel('To')
    colormap(flipud(colormap('gray')))
    
    title('Transitions')
    
    colorbar;
    
    state_count = sum(tmat,2);
    pmat = tmat./state_count;
    
    subplot(1,2,2)
    imagesc(pmat)
    xticks(1:7), xticklabels(states)
    yticks(1:7), yticklabels(states)
    ylabel('From')
    xlabel('To')
    colormap(flipud(colormap('gray')))
    % clim([0 1])
    colorbar
    title('Conditional Probability')

    rel_freq(:,i) = state_count(2:end)./sum(state_count(2:end));

%     suptitle(sessions{i})
end



% % This is not working
% beta_num = tmat./sum(tmat(:));
% beta_denom = prod(state_count/sum(state_count));
% B = log(beta_num ./ beta_denom);
% 
% subplot(1,2,2)
% imagesc(B)
% xticks(1:7), xticklabels(states)
% yticks(1:7), yticklabels(states)
% ylabel('From')
% xlabel('To')
% colormap(flipud(colormap('gray')))
% % clim([0 1])
% colorbar
% title('Transition Bias')

% figure
% p = pie(state_count(2:end), states(2:end));
% fontsize(gcf, 12, "points")
% arrayfun(@(data) set(data, 'FaceAlpha', 0.5), p(1:2:end))



% sort rel_freq so unilateral
%%
figure, 
b = bar(rel_freq);
b(1).FaceColor = 'k';
b(2).FaceColor = [0.3010 0.7450 0.9330];
xticklabels(states(2:end))
ylabel('Relative Frequency')
legend({'Spontaneous', 'Evoked'}, 'Location', 'NorthEast', 'Box', 'Off')
set(gca, 'FontSize', 12)



%% find multiple events within snippets

event_duration_all = cell(1,7);
for i = 1:length(event_duration)
    if ~isempty(event_duration{i})
        for j = 1:length(event_duration{i})
            event_duration_all{j} = cat(1,event_duration_all{j}, event_duration{i}{j});
        end
    else
        disp('hello')
    end
end


figure, hold on
for i = 1:length(event_duration_all)
    subplot(4, 2, i)
    if any(event_duration_all{i}<0)
        event_duration_all{i} = event_duration_all{i}(event_duration_all{i}>0);
    end
    [f,xi] = ksdensity(event_duration_all{i}); 
    plot(xi./fs, f)
    findpeaks(f, fs, 'MinPeakProminence', 0.005)
%     histogram(event_duration_all{i}./fs, 50);
%     vline(median(event_duration_all{i})/fs, 'r')
%     vline(2*median(event_duration_all{i})/fs, 'r:')
%     vline(3*median(event_duration_all{i})/fs, 'r:')
%     vline(4*median(event_duration_all{i})/fs, 'r:')
    xlim([0 1])
    xlabel('Duration (s)')
    title(labels(i))
    ylabel('Counts')
end



%%
% check if events within grooming episodes follows a pattern
longest_event = max(diff(idx));
b = nan(size(idx,2), longest_event, 3);
for ii = 1:size(b,1)
    tmp = bmat(idx(1, ii): idx(2, ii));
    ellip = tmp == 1;
    unilat = tmp == 2 | tmp == 3;
    bilat = tmp == 4;

    b(ii, 1:length(ellip), 1) = ellip;
    b(ii, 1:length(unilat), 2) = unilat;
    b(ii, 1:length(bilat), 3) = bilat;
end





%% network visualization

G = digraph(pmat,states);
figure, %plot(G)
p=plot(G, 'LineWidth', G.Edges.Weight*10);
numlinks = state_count*1;
p.MarkerSize = numlinks;
p.NodeFontSize = 12;
p.ArrowSize = 20;
p.EdgeAlpha = 0.25;
%%

%         for ii = 1:size(idx,1)
%             tmp = bmat2(idx(ii, 1): idx(ii, 2));
%             for jj = 1:length(tmp)
%                 if jj == 1 && tmp(jj) ~= 0
%                     tmat(1, tmp(jj)) = tmat(1, tmp(jj))+1;
%                     last_event = tmp(jj);
%                 else
%                      if tmp(jj)~=tmp(jj-1) && tmp(jj) ~= 1
%                          tmat(last_event, tmp(jj)) = tmat(last_event, tmp(jj)) + 1;
%                          last_event = tmp(jj);
%                      end
%                 end
%             end
%             tmat(last_event, 1) = tmat(last_event, 1) + 1;
%         
%         end
    
    
    
    %     figure,   plot(xt(bmat, fs), bmat), yticks(1:7), yticklabels(labels)
    % 
    %     bmat2 = bmat;
    %     bmat2(bmat2==3)=2;
    %     bmat2(bmat2==4)=3;
    %     bmat2(bmat2>4) = 0;
    %     figure,   plot(xt(bmat, fs), bmat2), yticks(1:3), yticklabels(["elliptical", "unilateral", "bilateral"])
    
    





%%

function [bmat2, idx] = aggregate(bmat, W, fs)
% aggregate grooming events by merging all events that occur within W
% seconds of each other. Then return list of start and stop indexes


if nargin < 2 || isempty(W), W = 3; end
if nargin < 3 || isempty(fs), fs = 90; end

W = round(W*fs); % convert seconds to samples

bmat2 = bmat>0;
bmat2 = movmax(bmat2, W);
bmat2 = movmin(bmat2, W);

% figure, plot(bmat>0), hold on, plot(bmat2-1)

time_on = find(diff(bmat2)==1)+1;
time_off = find(diff(bmat2)==-1)+1;

% edge cases where grooming occurs at the beginning or end
if bmat2(end), time_off = [time_off, length(bmat2)]; end
if bmat2(1), time_on = [1, time_on]; end

idx = [time_on; time_off]';
end


%%


clc, clear, %close all
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


current_mouse = '';

%%
k = 0.2;
for j = ger2_idx
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
                mvtOn{i} = mvtOn{i}(1:min_trial_length);
                lick_start{i} = lick_start{i}(1:min_trial_length);
                lick_timer{i} = lick_timer{i}(1:min_trial_length);
            end
            Vmaster = catcell(2, Vmaster);
            mvtOn = catcell(1, mvtOn)';
            lick_start = catcell(1, lick_start);
            lick_timer = catcell(1, lick_timer);
            new_trial = repmat([1 zeros(1, min_trial_length-1)], 1, num_trials);

            % ridge here
            disp('Building design matrix')
            opts.frameRate = fs;
            opts.sPostTime=round(fs*2);
            opts.mPreTime = ceil(0.5 * fs);  % precede motor events to capture preparatory activity in frames (used for eventType 3)
            opts.mPostTime = ceil(1 * fs);   % follow motor events for mPostStim in frames (used for eventType 3)
            opts.framesPerTrial = min_trial_length; % nr. of frames per trial
            opts.folds = 1; %nr of folds for cross-validation
            
            regressor_mat = [mvtOn' lick_start]; 
            % Full-Trial events:    new_trial
            % Peri-Stimulus events: 
            [dMat, regIdx] = makeDesignMatrix(regressor_mat, [3, 3], opts);
            regLabels = {'Movement', 'Lick', 'Timer'}; %some movement variables
%             [dMat, regIdx] = makeDesignMatrix(regressor_mat, [3, 3, 1], opts);
%             regLabels = {'Movment', 'Lick', 'Trial'}; %some movement variables
            fullR = [dMat, lick_timer];
            regIdx = [regIdx; max(regIdx)+1];

            disp('Running ridge regression with 10-fold cross-validation')
            [Vfull, fullBeta, ~, fullIdx, fullRidge, fullLabels] = crossValModel(fullR, Vmaster, regLabels, regIdx, regLabels, opts.folds);
%             save([fPath 'cvFull.mat'], 'Vfull', 'fullBeta', 'fullR', 'fullIdx', 'fullRidge', 'fullLabels', '-v7.3'); %save some results
            
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
            subplot(4, 1, 1)
            imagesc(reshape(fullMat, [128 128]))
            xticks([])
            c=colorbar;
            title('FullMat')
            c.Label.String = 'cvR^2';
            yticks([])
            for i = 1:length(regLabels)
                subplot(4, 1, i+1)
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
        boris_file = [data_dir, filesep, getAllFiles(data_dir, 'events.tsv')];
        dlc_speed_file = [data_dir, filesep, get_file_with_str(data_dir, 'speed.csv')];
        timestamp_file = [data_dir, filesep, get_file_with_str(data_dir, 'trim.txt')];
        ini_file = [data_dir, filesep, getAllFiles(data_dir, '.ini')];

        if ~any(isfile(brain_file) & isfile(boris_file) & isfile(dlc_speed_file) )
            continue
        end
    catch
        disp('Missing one or more of the required files. Skipping...')
        continue
    end
    % read ini file
    ini = ini2struct(ini_file);
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

    % Take off 2 frames due to possibility of dark frame at end creating
    % filter artifacts
    Vbrain = recastV(Umaster, U, s, V(:,1:end-2));
%     Vmaster{count} = Vbrain;
    [b, a] = butter(1, [0.01 10]/(fs/2));
    Vmaster{count} = filtfilt(b, a, Vbrain')';
    

    %% load DLC tracks
    disp('Loading DLC tracks')
    dlc_speed = readmatrix(dlc_speed_file);

    % consolidate limb speeds from all angles
    fl_l = sqrt(dlc_speed(:,1).^2 + dlc_speed(:,2).^2);
    fl_r = sqrt(dlc_speed(:,3).^2 + dlc_speed(:,4).^2);
    hl_l = dlc_speed(:,5);
    hl_r = dlc_speed(:,6);
    snout = sqrt(dlc_speed(:,7).^2 + dlc_speed(:,8).^2 + dlc_speed(:,9).^2);
    tailbase = sqrt(dlc_speed(:,10).^2 + dlc_speed(:,11).^2);
    
    all_movement = sqrt(fl_l.^2 + fl_r.^2 + hl_l.^2 + hl_r.^2 + snout.^2 + tailbase.^2);
    mvt = resample(all_movement, size(Vmaster{count},2), length(all_movement));
%     k = 10;
    mvtOn{count} = get_on_time(mvt > mean(mvt) + std(mvt));

    disp('Reading BORIS')
    [events, b_idx, b_tab] = read_boris(boris_file);
    lick = zeros(size(mvtOn{count}));
    lick(b_idx{1}(:,1)) = 1;
%     lick_start{count}(b_idx{1}(:,1)) = 1;
    lick_start{count} = lick;
    [~, lick_timer{count}] = start_timer(lick, k, fs);

%     lick_duration = b_idx{1}(end,end)-b_idx{1}(1,1);
%     lick_timer{count} = zeros(size(mvtOn{count}));
%     lick_timer{count}(b_idx{1}(1,1):b_idx{1}(end,end)) = 1:lick_duration+1;

    count = count + 1;

    clear U s V Vbrain
end


%%

fullLabels
visual = true;
cBetaRight = check_beta('Timer', fullLabels, fullIdx, Umaster, fullBeta{1}, Vfull, [], visual);

%%

test = reshape(Umaster*Vmaster, 128, 128, []);
test2 = reshape(Umaster*Vfull, 128, 128, []);
figure, plot(xt(test,fs), squeeze(mean(test, [1, 2])),'k'),
hold on
plot(xt(test, fs), squeeze(mean(test2, [1, 2]))), hold on,
% plot(xt(test,fs), lick_timer*.057)
%%

% figure

for i = 14
    
    test = reshape(Umaster*fullBeta{1}(fullIdx==i,:)', 128, 128, []);
    indices = floor(1:fs/4:size(test,3));
    if size(test,3)>1
        figure, imagesc(imtile(test, 'Frames', indices, 'GridSize', [1 length(indices)]))
        yticks([]);
        xticks((128:256:11*256)/2)
        xticklabels(-0.5:0.25:2)
        colormap(bluewhitered())
        c=colorbar;
        c.Label.String = 'Beta Kernel';
        xlabel('Time (s)')
    else
        figure, imagesc(test)
        xticks([])
        yticks([])
        caxis([-max(abs(clim)) max(abs(clim))]);
        colormap(bluewhitered())
        c=colorbar;
        c.Label.String = 'Beta';
        c.Label.FontSize = 14;
    end
end


%%


figure, plot(xt(lick_start, fs), lick_start*50, 'm')
hold on,
test = reshape(Umaster*Vmaster, 128, 128, []);
plot(xt(lick_start, fs), squeeze(mean(test, [1 2])),'k', 'LineWidth', 1)
xlabel('Time (s)')
axis tight
%%
%  start_timer(lick, k, fs)
test = aggregate(lick, k, fs)
%%

fs = 30;
k = 0.2;
[start_idx, b_timer] = start_timer(test, k, fs);

%%

function new_dat = get_on_time(data)
new_dat = zeros(size(data));
d = diff(data);
t_on = find(d>0)+ 1;

new_dat(t_on) = 1;

end

function [start_idx, behavior_timer] = start_timer(data, k, fs)
% takes a binary event vector (data) and aggregation window size (k)
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

%%

clc
% aa = catcell(1, grooming_regressors);
% bb = catcell(1, lick_regressor);

t = xt(grooming_regressors, fs);

cols = colororder('default');
figure, hold on
for i = 1:size(grooming_regressors,2)
    if ~any(grooming_regressors(:,i)), continue; end
%     plot(t, grooming_regressors(:,i)-i)
    patchplot(t(arr2idx(grooming_regressors(:,i))), [i-1 i], cols(i,:), 0.5)
end
% plot(t, lick_regressor - size(grooming_regressors,2)-1)
% plot(t, (activity_timer./max(activity_timer)) - size(grooming_regressors,2)-2)
% yticks(fliplr([0.5:1:size(grooming_regressors,2)+2]*-1))
% yticklabels(fliplr([grooming_behaviors, 'lick', 'timer']))



%%


clear, clc, close all
addpath('C:\Users\user\Documents\Nick\grooming\utils')
data_root = 'Y:\nick\behavior\grooming\1p';
mice = {'ECR2_thy1', 'GER2_ai94', 'HYL3_tTA', 'IBL2_tTA'};

for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    dff_path = [data_root, filesep, mice{j}, filesep, 'outputs'];
    dff_fig = getAllFiles(dff_path, 'dFF.fig');
%     dff_fig = 'summary.fig';
    % If there are multiple versions, use the most recent one
    if size(dff_fig,1) > 1
        dff_fig = sort(dff_fig);
        dff_fig = dff_fig{end};
    end
    h = openfig([dff_path, filesep, dff_fig]);
    for i = 1:length(h.Children)
        disp(h.Children(i).Title.String)
        switch(h.Children(i).Title.String)
            case 'RightMove' 
                rightmove(:,:,j) = h.Children(i).Children(end).CData;
            case 'LeftMove' 
                leftmove(:,:,j) = h.Children(i).Children(end).CData;
            case 'Large Bilateral'
                bilateral(:,:,j) = h.Children(i).Children(end).CData;
            case 'Elliptical'
                elliptical(:,:,j) = h.Children(i).Children(end).CData;
            case 'Elliptical Asymmetric'
                ellip_asymm(:,:,j) = h.Children(i).Children(end).CData;
            case 'Left Asymmetric'
                largeleft(:,:,j) = h.Children(i).Children(end).CData;
            case 'Right Asymmetric'
                largeright(:,:,j) = h.Children(i).Children(end).CData;
            case 'Left'
                left(:,:,j) = h.Children(i).Children(end).CData;
            case 'Right'
                right(:,:,j) = h.Children(i).Children(end).CData;
            case 'Lick'
                lick(:,:,j) = h.Children(i).Children(end).CData;
            otherwise
                continue
        end
    end

    close(h)
end


%%
% 
% figure
% for i = 1:4
% subplot(3,4,i)
% imagesc(leftlick(:,:,i))
% colorbar
% subplot(3,4,i+4), imagesc(lick(:,:,i));
% colorbar
% subplot(3,4,i+8), imagesc(left(:,:,i));
% colorbar
% end

%% overlay contours from diff mice
clc
% nanmask = zeros(128, 128, length(mice));

load('allenDorsalMap.mat');
clear nanmask
for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    atlas_tform = load([data_root, filesep, mice{j}, filesep, 'atlas_tform.mat']);
    warpmask = imwarp(mask, atlas_tform.tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
    nanmask(:,:,j) = warpmask;
%     nanmask(:,:,j) = mask;
end
nanmask(nanmask==0) = nan;





%%
clc, clear a v
tmp = [];
load('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat');

thresh = 85;
figure, axis off, hold on
for j = 1:length(mice)
    atlas_tform = load([data_root, filesep, mice{j}, filesep, 'atlas_tform.mat']);

    vars = ["lick", "right", "left", "elliptical", "largeright", "largeleft", "ellip_asymm"];%, "bilateral"];
         
    for i = 1:length(vars)    
        test = eval(vars(i));
        test = test(:,:,j) .* nanmask(:,:,j);
%         test = test.*nanmask(:,:,j);
%         test = imwarp(test, atlas_tform.tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));


        v = prctile(test(:), thresh);
        tmp = [tmp v];
%         disp(v)
        a{i}(:,:,j) = test >= v;

        subplot(1,length(vars), i), axis off, hold on
% subplot(1,1,1)
        for p = 1:length(dorsalMaps.edgeOutline)
            plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k');
        end
            if v > 0
                contourf(test, [v v], 'FaceAlpha', 0.25)
    
                title(vars(i));
            else
                disp('v is not greater than 0')
                disp(vars(i))
            end
            set(gca, 'YDir', 'reverse');
%         end
    end
%     legend(labs, 'Location', 'Best')
end

%% compute pairwise dice
for i = 1:length(vars)
    trialcount = 1;
    for j = 1:length(mice)
        for k = 1:length(mice)
            if k > j
                similarity(trialcount, i) = dice(a{i}(:,:,j), a{i}(:,:,k));
                trialcount = trialcount + 1;
            end
        end
    end
end


[p,tbl,stats] = anova1(similarity);
[c,m,h,gnames] = multcompare(stats);

figure, boxplot(similarity, 'Colors', 'k'),
xticklabels(vars)
ylabel('Pairwise Dice Similarity Coefficient')
ax = gca;
ax.FontSize = 12;

title(['One-Way ANOVA, p = ', num2str(p, 2)])




%% single trial dice matrix

clear, clc


fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = '';
trial = [];
mousecount = 0;

left = nan(128, 128, length(data_list{1}));
right = nan(128, 128, length(data_list{1}));
elliptical = nan(128, 128, length(data_list{1}));
largeleft = nan(128, 128, length(data_list{1}));
largeright = nan(128, 128, length(data_list{1}));
largebilateral = nan(128, 128, length(data_list{1}));
lick = nan(128, 128, length(data_list{1}));

for j = 1:length(data_list{1})
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
    if new_mouse       
        load([mouse_root_dir filesep 'mask.mat'])
        mask = double(mask);
        mask(mask==0) = nan;
        trialcount = 1;
        mousecount = [mousecount mousecount(end)+1];
    else
        trialcount = trialcount + 1;
        mousecount = [mousecount mousecount(end)];
    end
    trial = [trial, trialcount];
    
    dff_path = [data_dir filesep 'outputs'];
    h = openfig([dff_path, filesep, getAllFiles(dff_path, '_dFF.fig')]);
    current_mouse = mouse_id;

    vars = {'left', 'right', 'lick', 'elliptical', 'largeleft', 'largeright', 'largebilateral'};

    for i = 1:length(vars)
        for jj = 1:length(h.Children)
            if strcmp(h.Children(jj).Title.String, vars{i})
                switch vars{i}
                    case 'left'
                        left(:,:,j) = h.Children(jj).Children.CData.*mask;
                    case 'right'
                        right(:,:,j) = h.Children(jj).Children.CData.*mask;
                    case 'lick'
                        lick(:,:,j) = h.Children(jj).Children.CData.*mask;
                    case 'elliptical'
                        elliptical(:,:,j) = h.Children(jj).Children.CData.*mask;
                    case 'largeleft'
                        largeleft(:,:,j) = h.Children(jj).Children.CData.*mask;
                    case 'largeright'
                        largeright(:,:,j) = h.Children(jj).Children.CData.*mask;
                    case 'largebilateral'
                        largebilateral(:,:,j) = h.Children(jj).Children.CData.*mask;
                end
            end
        end
    end


    close(h)

end
mousecount = mousecount(2:end);

%%


combined_bstack = cat(3, lick, left, right, elliptical, largeleft, largeright, largebilateral);
for i = 1:size(combined_bstack, 3)
    tmp = combined_bstack(:,:,i);
    thresh = prctile(tmp(:), 90);
    combined_bstack(:,:,i) = tmp > thresh;
end

%%

for i = 1:size(combined_bstack,3)
    if sum(combined_bstack(:,:,i), [1 2]) == 0, continue; end
    for j = 1:size(combined_bstack,3)   
        if sum(combined_bstack(:,:,j), [1 2]) == 0, continue; end
        dice_matrix(i,j) = dice(combined_bstack(:,:,i), combined_bstack(:,:,j));
    end
end


%%



figure, imagesc(dice_matrix)
colormap(bluewhitered)
yticks(25:25:175)
xticks(25:25:175)


%%


data_root = 'Y:\nick\behavior\grooming\1p';
mice = {'ECR2_thy1', 'GER2_ai94', 'HYL3_tTA', 'IBL2_tTA'};

for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    dff_path = [data_root, filesep, mice{j}, filesep, 'outputs'];
    h = openfig([dff_path, filesep, getAllFiles(dff_path, '_dFF.fig')]);
    rightmove(:,:,j) = h.Children(8).Children.CData;
    leftmove(:,:,j) = h.Children(10).Children.CData;
    left(:,:,j) = h.Children(16).Children.CData;
    right(:,:,j) = h.Children(14).Children.CData;
    lick(:,:,j) = h.Children(12).Children.CData;
    elliptical(:,:,j) = h.Children(24).Children.CData;
    largeleft(:,:,j) = h.Children(22).Children.CData;
    largeright(:,:,j) = h.Children(20).Children.CData;
    dropright(:,:,j) = h.Children(2).Children.CData;
    bilateral(:,:,j) = h.Children(18).Children.CData;
    close(h)
end

%% overlay contours from diff mice
clc
nanmask = zeros(128, 128, length(mice));

for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    nanmask(:,:,j) = mask;
end
nanmask(nanmask==0) = nan;





%%
clc, clear a v
tmp = [];

thresh = 80;
figure, axis off, hold on
for j = 1:length(mice)

    vars = ["lick", "right", "left", "elliptical", "largeright", "largeleft", "bilateral"];
         
    for i = 1:length(vars)    
        test = eval(vars(i));
        test = test .* nanmask;
        test = test(:,:,j);

        v = prctile(test(:), thresh);
        tmp = [tmp v];
        disp(v)
        a{i}(:,:,j) = test >= v;

        subplot(1,length(vars), i), axis off, hold on
            if v > 0
                contourf(flipud(test), [v v], 'FaceAlpha', 0.25)
    
                title(vars(i));
            end
%         end
    end
%     legend(labs, 'Location', 'Best')
end

%% compute pairwise dice
for i = 1:length(vars)
    trialcount = 1;
    for j = 1:length(mice)
        for k = 1:length(mice)
            if k > j
                similarity(trialcount, i) = dice(a{i}(:,:,j), a{i}(:,:,k));
                trialcount = trialcount + 1;
            end
        end
    end
end




[p,tbl,stats] = anova1(similarity);
[c,m,h,gnames] = multcompare(stats);

figure, boxplot(similarity, 'Colors', 'k'),
xticklabels(vars)
ylabel('Pairwise Dice Similarity Coefficient')
ax = gca;
ax.FontSize = 12;

title(['One-Way ANOVA, p = ', num2str(p, 2)])




%% not using


%%   do the same thing on ridge residuals

clear, clc


fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
data_list = textscan(fileID, formatSpec);

    left_counter = 1;
    right_counter = 1;
    elliptical_counter = 1;
    lick_counter = 1;
    largeleft_counter = 1;
    largeright_counter = 1;
    largebilateral_counter = 1;


for j = 1:length(data_list{1})%1:length(data_list{1})+1
    % load experiment specific data into cell array
    ridge_dir = [data_list{1}{j}, filesep, 'ridge_outputs_video'];
    if ~isfolder(ridge_dir)
        disp('Ridge on Videp not performed yet for current experiment date')
        continue
    end
    
    h = openfig([ridge_dir, filesep, getAllFiles(ridge_dir, 'residuals')]);

    vars = {'left', 'right', 'lick', 'elliptical', 'largeleft', 'largeright', 'largebilateral'};

    for i = 1:length(vars)
        for jj = 1:length(h.Children)
            if strcmp(h.Children(jj).Title.String, vars{i})
                switch vars{i}
                    case 'left'
                        left(:,:,left_counter) = h.Children(jj).Children.CData;
                        left_counter = left_counter + 1;
                    case 'right'
                        right(:,:,right_counter) = h.Children(jj).Children.CData;
                        right_counter = right_counter + 1;
                    case 'lick'
                        lick(:,:,lick_counter) = h.Children(jj).Children.CData;
                        lick_counter = lick_counter + 1;
                    case 'elliptical'
                        elliptical(:,:,elliptical_counter) = h.Children(jj).Children.CData;
                        elliptical_counter = elliptical_counter + 1;
                    case 'largeleft'
                        largeleft(:,:,largeleft_counter) = h.Children(jj).Children.CData;
                        largeleft_counter = largeleft_counter + 1;
                    case 'largeright'
                        largeright(:,:,largeright_counter) = h.Children(jj).Children.CData;
                        largeright_counter = largeright_counter + 1;
                    case 'largebilateral'
                        largebilateral(:,:,largebilateral_counter) = h.Children(jj).Children.CData;
                        largebilateral_counter = largebilateral_counter + 1;
                end
            end
        end
    end

    close(h)

end


%%

thresh = 95;
figure, hold on
for j = 1:length(vars)
    subplot(1, length(vars), j), hold on
    tmp0 = eval(vars{j});
%     tmp = mean(tmp0, 3);
    v = prctile(tmp(:), thresh);
%     imagesc(flipud(tmp)), colorbar

%     contourf(flipud(tmp), [v v], 'FaceAlpha', 0.25)
    for i = 1:size(tmp0,3)
%         tmp = abs(tmp0(:,:,i));
        tmp = -tmp0(:,:,i);
        v = prctile(tmp(:), thresh);
        contourf(flipud(tmp), [v v], 'FaceAlpha', 0.25)
        
    end
    title(vars{j})
end

%%
clear, clc

data_root = 'Y:\nick\behavior\grooming\1p';
mice = {'ECR2_thy1', 'GER2_ai94', 'HYL3_tTA', 'IBL2_tTA'};



for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    dff_path = [data_root, filesep, mice{j}, filesep, 'outputs'];
    h = openfig([dff_path, filesep, getAllFiles(dff_path, 'ridge_summary.fig')]);
    rightmove(:,:,j) = h.Children(2).Children.CData;
    leftmove(:,:,j) = h.Children(4).Children.CData;
    left(:,:,j) = h.Children(10).Children.CData;
    right(:,:,j) = h.Children(8).Children.CData;
    lick(:,:,j) = h.Children(6).Children.CData;
    elliptical(:,:,j) = h.Children(18).Children.CData;
    largeleft(:,:,j) = h.Children(16).Children.CData;
    largeright(:,:,j) = h.Children(14).Children.CData;
    bilateral(:,:,j) = h.Children(12).Children.CData;
    close(h)
end

%% overlay contours from diff mice
clc
nanmask = zeros(128, 128, length(mice));

for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    nanmask(:,:,j) = mask;
end
nanmask(nanmask==0) = nan;





%%
clc, clear a v
tmp = [];

thresh = 80;
figure, axis off, hold on
for j = 1:length(mice)

    vars = ["lick", "right", "left", "elliptical", "largeright", "largeleft", "bilateral"];
         
    for i = 1:length(vars)    
        test = eval(vars(i));
        test = test .* nanmask;
        test = test(:,:,j);

        v = prctile(test(:), thresh);
        a{i}(:,:,j) = test >= v;

        subplot(1,length(vars), i), axis off, hold on
            if v > 0
                contourf(flipud(test), [v v], 'FaceAlpha', 0.25)
    
                title(vars(i));
            end
%         end
    end
%     legend(labs, 'Location', 'Best')
end

%% compute pairwise dice
for i = 1:length(vars)
    trialcount = 1;
    for j = 1:length(mice)
        for k = 1:length(mice)
            if k > j
                similarity(trialcount, i) = dice(a{i}(:,:,j), a{i}(:,:,k));
                trialcount = trialcount + 1;
            end
        end
    end
end




[p,tbl,stats] = anova1(similarity);
[c,m,h,gnames] = multcompare(stats);

figure, boxplot(similarity, 'Colors', 'k'),
xticklabels(vars)
ylabel('Pairwise Dice Similarity Coefficient')
ax = gca;
ax.FontSize = 12;

title(['One-Way ANOVA, p = ', num2str(p, 2)])






%%   do the same thing on ridge unique explained var

clear, clc

data_root = 'Y:\nick\behavior\grooming\1p';
mice = {'ECR2_thy1', 'GER2_ai94', 'HYL3_tTA', 'IBL2_tTA'};



for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    dff_path = [data_root, filesep, mice{j}, filesep, 'outputs'];
    h = openfig([dff_path, filesep, getAllFiles(dff_path, 'ridge_summary.fig')]);
    rightmove(:,:,j) = h.Children(2).Children.CData;
    leftmove(:,:,j) = h.Children(4).Children.CData;
    left(:,:,j) = h.Children(10).Children.CData;
    right(:,:,j) = h.Children(8).Children.CData;
    lick(:,:,j) = h.Children(6).Children.CData;
    elliptical(:,:,j) = h.Children(18).Children.CData;
    largeleft(:,:,j) = h.Children(16).Children.CData;
    largeright(:,:,j) = h.Children(14).Children.CData;
    bilateral(:,:,j) = h.Children(12).Children.CData;
    close(h)
end

%% overlay contours from diff mice
clc
nanmask = zeros(128, 128, length(mice));

for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    nanmask(:,:,j) = mask;
end
nanmask(nanmask==0) = nan;





%%
clc, clear a v
tmp = [];

thresh = 80;
figure, axis off, hold on
for j = 1:length(mice)

    vars = ["lick", "right", "left", "elliptical", "largeright", "largeleft", "bilateral"];
         
    for i = 1:length(vars)    
        test = eval(vars(i));
        test = test .* nanmask;
        test = test(:,:,j);

        v = prctile(test(:), thresh);
        a{i}(:,:,j) = test >= v;

        subplot(1,length(vars), i), axis off, hold on
            if v > 0
                contourf(flipud(test), [v v], 'FaceAlpha', 0.25)
    
                title(vars(i));
            end
%         end
    end
%     legend(labs, 'Location', 'Best')
end

%% compute pairwise dice
for i = 1:length(vars)
    trialcount = 1;
    for j = 1:length(mice)
        for k = 1:length(mice)
            if k > j
                similarity(trialcount, i) = dice(a{i}(:,:,j), a{i}(:,:,k));
                trialcount = trialcount + 1;
            end
        end
    end
end




[p,tbl,stats] = anova1(similarity);
[c,m,h,gnames] = multcompare(stats);

figure, boxplot(similarity, 'Colors', 'k'),
xticklabels(vars)
ylabel('Pairwise Dice Similarity Coefficient')
ax = gca;
ax.FontSize = 12;

title(['One-Way ANOVA, p = ', num2str(p, 2)])


