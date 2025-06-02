clc, clear, %close all
fileID = fopen('expt3_datalist.txt','r');
% addpath('C:\Users\user\Documents\Nick\grooming\utils')

% addpath('C:\Users\user\Documents\Nick\ridgeModel');
% addpath('C:\Users\user\Documents\Nick\ridgeModel\widefield')
% addpath('C:\Users\user\Documents\Nick\ridgeModel\smallStuff') 
% addpath('C:\Users\user\Documents\Nick\grooming\utils')
addpath('/home/user/Documents/grooming/utils')
addpath('/home/user/Documents/grooming/ridgeModel')
addpath('/home/user/Documents/grooming/ridgeModel/widefield')
addpath('/home/user/Documents/grooming/ridgeModel/smallStuff')
    
formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
data_list = data_list{1};

ger2_idx = find(contains(data_list, 'GER2'));
hyl3_idx = find(contains(data_list, 'HYL3'));
ecr2_idx = find(contains(data_list, 'ECR2'));
gt33_idx = find(contains(data_list, 'GT33'));

current_mouse = '';

%%
k = 0.2;
expts_to_analyze = 1:length(data_list);
for j = 1:length(expts_to_analyze)+1
     try
        data_dir = data_list{expts_to_analyze(j)};
        if isunix % working on linux computer - modify paths
            data_dir = strrep(data_dir, '\', '/');
            data_dir = strrep(data_dir, 'Y:/', '/media/user/teamshare/TM_Lab/');
            % data_dir = strrep(data_dir, 'Y:/', '/home/user/teamshare/TM_Lab/');
        end  
        disp(['Starting ' data_dir])
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
                mvt{i} = mvt{i}(1:min_trial_length);
                mvtOn{i} = mvtOn{i}(1:min_trial_length);
                lick_start{i} = lick_start{i}(1:min_trial_length);
                lick_rate{i} = lick_rate{i}(1:min_trial_length);
            end
            Vmaster = catcell(2, Vmaster);
            mvt = catcell(1, mvt);
            mvtOn = catcell(1, mvtOn)';
            lick_start = catcell(1, lick_start);
            lick_rate = catcell(1, lick_rate);


% commented for speed
            % ridge here
            disp('Building design matrix')
            opts.frameRate = fs;
            opts.sPostTime=round(fs*2);
            opts.mPreTime = ceil(0.5 * fs);  % precede motor events to capture preparatory activity in frames (used for eventType 3)
            opts.mPostTime = ceil(2 * fs);   % follow motor events for mPostStim in frames (used for eventType 3)
            opts.framesPerTrial = min_trial_length; % nr. of frames per trial
            opts.folds = 10; %nr of folds for cross-validation
            
            regressor_mat = [mvtOn' lick_start]; 
            % Full-Trial events:    new_trial
            % Peri-Stimulus events: 
            [dMat, regIdx] = makeDesignMatrix(regressor_mat, [3, 3], opts);
            regLabels = {'Movement', 'Lick'}; %some movement variables
            fullR = [dMat];

            
            disp('Running ridge regression with 10-fold cross-validation')
            [Vfull, fullBeta, ~, fullIdx, fullRidge, fullLabels] = crossValModel(fullR, Vmaster, regLabels, regIdx, regLabels, opts.folds);
            save([fPath 'cvFull.mat'], 'Vfull', 'fullBeta', 'fullR', 'fullIdx', 'fullRidge', 'fullLabels', '-v7.3'); %save some results
            
            fullMat = modelCorr(Vmaster,Vfull,Umaster) .^2;

            % break

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
            subplot(3, 1, 1)
            imagesc(reshape(fullMat, [128 128]))
            xticks([])
            c=colorbar;
            title('FullMat')
            c.Label.String = 'cvR^2';
            yticks([])
            for i = 1:length(regLabels)
                subplot(3, 1, i+1)
                imagesc(reshape(fullMat - reducedMat(:,i), [128 128]))
                c=colorbar;
                title(regLabels{i})
                c.Label.String = '\DeltaR^2';
            
                xticks([])
                yticks([])
            end

            savefig(gcf, [fPath 'drinking_summary.fig'])
% commented for speed ^^


            % make averaged map
            dff = zscore(reshape(Umaster*Vmaster, 128, 128, size(Vmaster,2)), [], 3);
            figure, imagesc(mean(dff(:,:,lick_rate>=1),3).*mask)
            colorbar
            colormap(bluewhitered())
            nanmask = mask;
            nanmask(mask==0) = nan;
            savefig(gcf, [fPath 'average_licking_map.fig'])

            save([fPath 'outputs.mat'], 'Vmaster', 'lick*', 'mvt*', 'min_trial_length', 'nanmask', '-v7.3')            

            % all mice completed - break the loop
            if j == length(data_list{1})+1, break; end
        end        
        count = 1;
        disp('Loading master basis set')
        load([mouse_root_dir filesep 'Umaster.mat'])
        clear Vmaster mvtOn lick_start lick_timer lick_rate mvt 
        current_mouse = mouse_id;
        load([mouse_root_dir filesep 'mask.mat'])

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
    mov = resample(all_movement, size(Vmaster{count},2), length(all_movement));
    mvt{count} = mov > mean(mov) + std(mov);
    mvtOn{count} = get_on_time(mvt{count});

    disp('Reading BORIS')
    [events, b_idx, b_tab] = read_boris(boris_file);
    lick = zeros(size(mvtOn{count}));
    lick(b_idx{1}(:,1)) = 1;
    lick_rate{count} = movsum(lick, fs);
    lick_start{count} = lick;


    count = count + 1;

    clear U s V Vbrain
end





%%

fullLabels
visual = true;
cBetaRight = check_beta('Lick', fullLabels, fullIdx, Umaster, mean(catcell(3,fullBeta),3), Vfull, [], visual);

%%

%%

load('allenDorsalMap.mat');
fPath = 'Y:\nick\behavior\grooming\figures';

data_list = {'Y:\pankaj\closedloop_rig5_data\GER2_ai94_drinking', ...
    'Y:\pankaj\closedloop_rig5_data\HYL3_tta_drinking', ...
    'Y:\pankaj\closedloop_rig5_data\ECR2_thy1_drinking', ...
    'Y:\pankaj\closedloop_rig5_data\GT33_tta_drinking'};

example_mouse = 2;
figure, hold on
clear lick_map
axis off, hold on
for p = 1:length(dorsalMaps.edgeOutline)
    plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k', 'LineWidth', 1);
    xticks([])
    yticks([])
end
set(gca, 'YDir', 'reverse');
for i = 1:length(data_list)
    load([data_list{i}, filesep, 'mask.mat'])
    h = openfig([data_list{i}, filesep, 'outputs' filesep, 'average_licking_map.fig']);
    load([data_list{i}, filesep, 'atlas_tform.mat'])

    nanmask = mask;
    nanmask(mask==0) = nan;
    lick_map(:,:,i) = imwarp(h.Children(2).Children.CData.*nanmask, tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)), 'FillValues', nan);

    if i == example_mouse
        img_data = h.Children(2).Children.CData.*nanmask;
        figure; bb=imagesc(img_data);
        set(bb, 'AlphaData', ~isnan(img_data))
        caxis([-0.2 1.2])
        axis off
        exportgraphics(gcf, fix_path([fPath, filesep, 'drinking_lick_dff_example.png']), 'Resolution', 300)
        close(gcf)
        figure, axis off, caxis([-0.2 1.2])
        colorbar
        exportgraphics(gcf, fix_path([fPath, filesep, 'drinking_lick_dff_example_colorbar.png']), 'Resolution', 300)
        close(gcf)
    end
    close(h)
    v = prctile(lick_map(:,:,i), 80);
    binary_maps(:,:,i) = lick_map(:,:,i)>=v;
    contourf(binary_maps(:,:,i).*i, [i-0.1 i-0.1], 'FaceAlpha', 0.25)

    
    
end

% exportgraphics(gcf, fix_path([fPath, filesep, 'drinking_lick_dff_contours.png']), 'Resolution', 300)
%% plot contours for lick dcvr2

figure, hold on
clear lick_map binary_maps
axis off, hold on
for p = 1:length(dorsalMaps.edgeOutline)
    plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k', 'LineWidth', 1);
    xticks([])
    yticks([])
end
set(gca, 'YDir', 'reverse');
for i = 1:length(data_list)
    load([data_list{i}, filesep, 'mask.mat'])
    nanmask = mask;
    nanmask(mask==0) = nan;
    load([data_list{i}, filesep, 'atlas_tform.mat'])
    h = openfig([data_list{i}, filesep, 'outputs' filesep, 'drinking_summary.fig']);    
    for j = 1:length(h.Children)
        switch(h.Children(j).Title.String)
            case 'Lick' 
                lick_map(:,:,i) = imwarp(h.Children(j).Children.CData.*nanmask, tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)), 'FillValues', nan);
            otherwise
                continue
        end
    end

    close(h)

    v = prctile(lick_map(:,:,i), 80);
    binary_maps(:,:,i) = lick_map(:,:,i)>=v;
    contourf(binary_maps(:,:,i).*i, [i-0.1 i-0.1], 'FaceAlpha', 0.25)

end

% exportgraphics(gcf, fix_path([fPath, filesep, 'drinking_lick_dcvr2_contours.png']), 'Resolution', 300)


%% plot contours for move dcvr2
fPath = 'Y:\nick\behavior\grooming\figures';
load('allenDorsalMap.mat');

figure, hold on
clear move_map binary_maps
axis off, hold on
for p = 1:length(dorsalMaps.edgeOutline)
    plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k', 'LineWidth', 1);
    xticks([])
    yticks([])
end
set(gca, 'YDir', 'reverse');
for i = 1:length(data_list)
    load([data_list{i}, filesep, 'mask.mat'])
    nanmask = mask;
    nanmask(mask==0) = nan;
    load([data_list{i}, filesep, 'atlas_tform.mat'])
    h = openfig([data_list{i}, filesep, 'outputs' filesep, 'drinking_summary.fig']);    
    for j = 1:length(h.Children)
        switch(h.Children(j).Title.String)
            case 'Movement' 
                move_map(:,:,i) = imwarp(h.Children(j).Children.CData.*nanmask, tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)), 'FillValues', nan);
            otherwise
                continue
        end
    end

    close(h)

    v = prctile(move_map(:,:,i), 80);
    binary_maps(:,:,i) = move_map(:,:,i)>=v;
    contourf(binary_maps(:,:,i).*i, [i-0.1 i-0.1], 'FaceAlpha', 0.25)

end

% exportgraphics(gcf, fix_path([fPath, filesep, 'drinking_move_dcvr2_contours.png']), 'Resolution', 300)


%% Plot example trials
% use first 4 from hyl3

clear, clc
data_list = {'Y:\pankaj\closedloop_rig5_data\GER2_ai94_drinking', ...
    'Y:\pankaj\closedloop_rig5_data\HYL3_tta_drinking', ...
    'Y:\pankaj\closedloop_rig5_data\ECR2_thy1_drinking', ...
    'Y:\pankaj\closedloop_rig5_data\GT33_tta_drinking'};

example_mouse = data_list{2};
load([example_mouse, filesep, 'outputs', filesep, 'outputs.mat'])
load([example_mouse, filesep, 'Umaster.mat'])
load([example_mouse, filesep, 'mask.mat'])
nanmask = mask;
nanmask(nanmask==0) = nan;

num_trials = 4;
figure, 
for i = 1:num_trials
    idx = min_trial_length*(i-1)+1:min_trial_length*i;
    dff = zscore(reshape(Umaster*Vmaster(:,idx), 128, 128, []),[],3);
    t = xt(dff, 30, 3);
    subplot(num_trials,1,i), hold on, 
    plot(t, squeeze(mean(dff.*nanmask, [1 2], 'omitnan')), 'k', 'LineWidth', 1)
    patchplot(t(arr2idx(lick_rate(idx)>=1)), [-4 4], 'm', 0.2)
    patchplot(t(arr2idx(mvt(idx))), [-4 4], 'b', 0.2)
    plot(t, squeeze(mean(dff.*nanmask, [1 2], 'omitnan')), 'k', 'LineWidth', 1)
%     plot(t, smoothdata(lick_rate(idx), 'gaussian', 10))
    axis([0 t(end) -4 4])
    line([0 t(end)], [0 0], 'Color', [0 0 0])
    
end

% saveas(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'drinking_4trials.svg']))
%% Plot lick rates for trial structure figure.
% Use just the example mouse since others have variable durations

for i = 1:12
    idx = min_trial_length*(i-1)+1:min_trial_length*i;
    lick_rate_compile(i,:) = lick_rate(idx);
end
% idx = 1:4*min_trial_length;
% test = zscore(reshape(Umaster*Vmaster(:,idx), 128, 128, []),[],3);
% figure
% imagesc(mean(test(:,:,lick_rate(idx)>=1),3))
figure, shadedErrorBar(t, mean(lick_rate_compile), std(lick_rate_compile)./sqrt(size(lick_rate_compile,1)))
axis tight
% figure, shadedErrorBar(t, mean(lick_rate_compile), std(lick_rate_compile))
xlabel('Time (s)')
ylabel('Licks per second')

%% Plot the averaged dFF map
figure
h = openfig('Y:\pankaj\closedloop_rig5_data\HYL3_tta_drinking\outputs\average_licking_map.fig')
colormap parula


%% Plot betas for Lick and movement

data_list = {'Y:\pankaj\closedloop_rig5_data\GER2_ai94_drinking', ...
    'Y:\pankaj\closedloop_rig5_data\HYL3_tta_drinking', ...
    'Y:\pankaj\closedloop_rig5_data\ECR2_thy1_drinking', ...
    'Y:\pankaj\closedloop_rig5_data\GT33_tta_drinking'};

example_mouse = data_list{3};

load([example_mouse, filesep, 'outputs', filesep, 'cvFull.mat'])
load([example_mouse, filesep, 'Umaster.mat'])
load([example_mouse, filesep, 'mask.mat'])
load([example_mouse filesep 'atlas_tform.mat'])
fPath = [example_mouse filesep 'outputs' filesep];
[seeds, labels] = get_seeds();
close()



nanmask = mask;
nanmask(nanmask==0) = nan;

%%
visual = false;
cBetaRight = check_beta('Lick', fullLabels, fullIdx, Umaster, mean(catcell(3,fullBeta),3), Vfull, [], visual);

beta = cBetaRight.*nanmask;
figure('Position', [45 593 602 180]), 
img_data = [beta(:,:,15), beta(:,:,30), beta(:,:,45)];
bb = imagesc(img_data);
set(bb, 'AlphaData', ~isnan(img_data))
caxis([-5 5])
axis off
hold on
rois = {'MOS_2-L', 'SSP-ul-L', 'PTLp-L'};

cmap = colororder();
for i = 1:length(rois)
    tmp = seeds(strcmp(labels, rois(i)),:);
    % seeds are in atlas coords. Use inverse of atlas tform to get them
    % back in original space

    [u(i), v(i)] = transformPointsInverse(tform, tmp(1), tmp(2));
    plot(u(i),v(i), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(i,:))
end

exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures',filesep, 'drinking_lick_betas_montage.png']), 'Resolution', 300)

figure, axis off, caxis([-5 5]), colorbar
exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures',filesep, 'drinking_licks_betas_montage_colorbar.png']), 'Resolution', 300)


figure, hold on
tmp = getTimeseries(beta, [u' v']);
t = xt(tmp, 30, 1)-0.5;
plot(t, tmp)
line([t(1) t(end)], [0 0], 'Color', [0 0 0])
vline([0 0.5 1], {'k', 'k', 'k'})
% saveas(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'drinking_betas_timeseries.svg']))

%% movement betas
visual = false;
cBetaRight = check_beta('Movement', fullLabels, fullIdx, Umaster, mean(catcell(3,fullBeta),3), Vfull, [], visual);


beta = cBetaRight.*nanmask;
figure('Position', [45 593 602 180]), 
img_data = [beta(:,:,15), beta(:,:,30), beta(:,:,45)];
bb = imagesc(img_data);
set(bb, 'AlphaData', ~isnan(img_data))
caxis([-5 20])

hold on
rois = {'MOS_2-L', 'SSP-ll-L', 'PTLp-L'};

cmap = colororder();
for i = 1:length(rois)
    tmp = seeds(strcmp(labels, rois(i)),:);
    % seeds are in atlas coords. Use inverse of atlas tform to get them
    % back in original space

    [u(i), v(i)] = transformPointsInverse(tform, tmp(1), tmp(2));
    plot(u(i),v(i), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(i,:))
end
axis off
exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures',filesep, 'drinking_movement_betas_montage.png']), 'Resolution', 300)

figure, axis off, caxis([-5 20]), colorbar
exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures',filesep, 'drinking_movement_betas_montage_colorbar.png']), 'Resolution', 300)

figure, hold on
tmp = getTimeseries(beta, [u' v']);
t = xt(tmp, 30, 1)-0.5;
plot(t, tmp)
line([t(1) t(end)], [0 0], 'Color', [0 0 0])
vline([0 0.5 1], {'k', 'k', 'k'})
% saveas(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'drinking_movement_betas_timeseries.svg']))



function new_dat = get_on_time(data)
new_dat = zeros(size(data));
d = diff(data);
t_on = find(d>0)+ 1;

new_dat(t_on) = 1;

end
