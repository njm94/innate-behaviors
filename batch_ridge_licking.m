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
            end
            Vmaster = catcell(2, Vmaster);
            mvtOn = catcell(1, mvtOn)';
            lick_start = catcell(1, lick_start);
            new_trial = repmat([1 zeros(1, min_trial_length-1)], 1, num_trials);

            % ridge here
            disp('Building design matrix')
            opts.frameRate = fs;
            opts.sPostTime=round(fs*2);
            opts.mPreTime = ceil(0.5 * fs);  % precede motor events to capture preparatory activity in frames (used for eventType 3)
            opts.mPostTime = ceil(2 * fs);   % follow motor events for mPostStim in frames (used for eventType 3)
            opts.framesPerTrial = min_trial_length; % nr. of frames per trial
            opts.folds = 1; %nr of folds for cross-validation
            
            regressor_mat = [mvtOn' lick_start]; 
            % Full-Trial events:    new_trial
            % Peri-Stimulus events: 
            [dMat, regIdx] = makeDesignMatrix(regressor_mat, [3, 3], opts);
            regLabels = {'Movment', 'Lick'}; %some movement variables
%             [dMat, regIdx] = makeDesignMatrix(regressor_mat, [3, 3, 1], opts);
%             regLabels = {'Movment', 'Lick', 'Trial'}; %some movement variables
            fullR = [dMat];

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
            subplot(3, 3, 1)
            imagesc(reshape(fullMat, [128 128]))
            xticks([])
            c=colorbar;
            title('FullMat')
            c.Label.String = 'cvR^2';
            yticks([])
            for i = 1:length(regLabels)
                subplot(3, 3, i+1)
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
    
%     try ini = ini.sentech_give_rewards;
%     catch 
%         ini = ini.sentech_dlc_live;
%     end
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
    % Vmaster{count} = Vbrain;
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
    k = 10;
    mvtOn{count} = get_on_time(movmax(mvt > mean(mvt) + std(mvt), k));

    disp('Reading BORIS')
    [events, b_idx, b_tab] = read_boris(boris_file);
    lick_start{count} = zeros(size(mvtOn{count}));
    lick_start{count}(b_idx{1}(:,1)) = 1;

    count = count + 1;

    clear U s V Vbrain
end


%%

clc
visual = true;
cBetaRight = check_beta('Lick', fullLabels, fullIdx, Umaster, fullBeta{1}, Vfull, [], visual);

%%

figure
indices = floor(1:7.5:76);
for i = 1:2
    
    test = reshape(Umaster*fullBeta{1}(fullIdx==i,:)', 128, 128, []);
    subplot(2,1,i), imagesc(imtile(test, 'Frames', indices, 'GridSize', [1 length(indices)]))
    yticks([]);
    xticks((128:256:11*256)/2)
    xticklabels(-0.5:0.25:2)
    colormap(bluewhitered())
    c=colorbar;
    c.Label.String = 'Beta Kernel';
    xlabel('Time (s)')
end

%%

function new_dat = get_on_time(data)
new_dat = zeros(size(data));
d = diff(data);
t_on = find(d>0)+ 1;

new_dat(t_on) = 1;

end