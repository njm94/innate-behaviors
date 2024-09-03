% Make behavior features for UMAP representation
%
% 1) Read boris behavior events
% 2) Read deeplabcut file
% 3) Create NxM feature matrix, where N is number of facial grooming
% strokes and M is number of features
% 4) Create Nx1 label matrix
% 5) Compute UMAP on behavior feature matrix color-coded by manual labels
%
% Features to use:
%   Duration (num frames)
%   Interspersed with lick (0 or 1)
%   Depth in chain (0 to n) ??
%   Left paw Ymax / Right paw Ymax (ratio)
%   Paw starting position relative to nose
%       Left,   x
%       Left,   y
%       Right,  x
%       Right,  y
%   




%%

clc, clear
fileID = fopen('expt1_datalist.txt','r');
addpath('C:\Users\user\Documents\Nick\grooming\utils')

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
mp_list = {'Y:\nick\behavior\grooming\2p\ETR2_thy1\20231113143925'; ...
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231113155903'; ...
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231115174148';
    };
mp_list = {'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240729'; ...
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240731';
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240802'};
data_list{1} = [data_list{1}; mp_list];
data_list{1} = mp_list;
current_mouse = '';

fs = 90 ;

thy1_idx = 1:7;
ai94_idx = 8:13;
camk_idx = 14:25;
mp_idx = 26:28;


%% 

% initialize feature and label matrices
fMat = []; 
lMat = [];


behavior_duration = [];
behavior_label = [];
l_paw_start_x = [];
l_paw_start_y = [];
r_paw_start_x = [];
r_paw_start_y = [];
l_paw_yrange = [];
l_paw_xrange = [];
r_paw_yrange = [];
r_paw_xrange = [];
l_paw_max_y = [];
r_paw_max_y = [];
lick_interspersed = [];
fll_distance = [];
flr_distance = [];
imm_events = [];
fll_max_speed = [];
flr_max_speed = [];
fll_avg_speed = [];
flr_avg_speed = [];

%%


for j = 1:length(data_list{1})%mp_idx %23:length(data_list{1})+1
     try
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
    catch
        new_mouse = true;
    end
    if new_mouse 
        if ~isempty(current_mouse)
            % if all the trials have completed for the previous mouse,
            % compile all behavior frames and take the average

            if j == length(data_list{1})+1, break; end
        end        
        count = 1;
        current_mouse = mouse_id;

        fPath = [mouse_root_dir filesep 'outputs' filesep];       
    end

    try
        dlc_pos_file = [data_dir, filesep, get_file_with_str(data_dir, '1030000.csv')];
        dlc_speed_file = [data_dir, filesep, get_file_with_str(data_dir, 'speed.csv')];
        boris_file = [data_dir, filesep, get_file_with_str(data_dir, 'events.tsv')];
    catch
        disp('Missing one or more of the required files. Skipping...')
        continue
    end

    % load grooming events from BORIS file
    [events, b_idx, boris] = read_boris(boris_file);

    % consolidate lick events
    lick_events = events(:,contains(events.Properties.VariableNames, 'Lick'));
    lick_events = any(table2array(lick_events),2);

    % remove lick and point events from event matrix
    idx = contains(events.Properties.VariableNames, 'Lick') | ...
        contains(events.Properties.VariableNames, 'Drop') | ...
        contains(events.Properties.VariableNames, 'Video');
    stroke_events = removevars(events, idx);

    % consolidate stroke_events to get aggregated grooming chains to assess
    % depth in chain as a feature
    consolidated_strokes = aggregate(any(table2array(stroke_events),2), 3);  


    %% load DLC tracks
    disp('Loading DLC tracks')
    dlc_pos = readmatrix(dlc_pos_file);
    dlc_pos = dlc_pos(1:size(stroke_events,1), :);

    nose_x = median(dlc_pos(:,1));
    nose_y = median(dlc_pos(:,2));
    
    flr_x = dlc_pos(:,4);
    flr_y = dlc_pos(:,5);
    
    fll_x = dlc_pos(:,7);
    fll_y = dlc_pos(:,8);

    dlc_speed = readmatrix(dlc_speed_file);
    fll_speed = dlc_speed(:,1);
    flr_speed = dlc_speed(:,2);

    %%
    figure
    for i = 1:size(stroke_events,2)
        tmp = logical(table2array(stroke_events(:,i)));
        subplot(1,7,i), hold on
        plot(flr_x(tmp), -flr_y(tmp))
        plot(fll_x(tmp), -fll_y(tmp))
        axis([200 450 -250 0])
        title(stroke_events.Properties.VariableNames{i})
    end

%     figure, plot(flr_x(table2array(stroke_events(:,1)), -flr_y(table2array(stroke_events(:,1)))))


%%

%   Duration (num frames)
%   Interspersed with lick (0 or 1)
%   Depth in chain (0 to n) ??
%   Left paw Ymax / Right paw Ymax (ratio)
%   Paw starting position relative to nose
%       Left,   x
%       Left,   y
%       Right,  x
%       Right,  y

    %% 

    stroke_types = stroke_events.Properties.VariableNames;

    for i = 1:length(stroke_types)
        idx = strcmpi(stroke_types{i}, boris.Behavior);
        behavior_frames = boris.ImageIndex(idx);

        
        start_idx = behavior_frames(1:2:end);
        stop_idx = behavior_frames(2:2:end);

        if length(start_idx) > 1
            % flag events that occur immediately one after the other
            frame_diff_between_events = start_idx(2:end)-stop_idx(1:end-1);
            flag_events = frame_diff_between_events < 3;
            if flag_events(end)
                flag_events = [flag_events; 1];
            else
                flag_events = [flag_events; 0];
            end
        else
            flag_events = 0;
        end

        imm_events = cat(1, imm_events, flag_events);


        l_paw_start_x = cat(1, l_paw_start_x, nose_x-fll_x(start_idx));
        l_paw_start_y = cat(1, l_paw_start_y, nose_y-fll_y(start_idx));
        r_paw_start_x = cat(1, r_paw_start_x, nose_x-flr_x(start_idx));
        r_paw_start_y = cat(1, r_paw_start_y, nose_y-flr_y(start_idx));

        for j = 1:length(start_idx)
            l_paw_yrange = cat(1, l_paw_yrange, range(fll_y(start_idx(j):stop_idx(j))));
            l_paw_xrange = cat(1, l_paw_xrange, range(fll_x(start_idx(j):stop_idx(j))));
            r_paw_yrange = cat(1, r_paw_yrange, range(flr_y(start_idx(j):stop_idx(j))));
            r_paw_xrange = cat(1, r_paw_xrange, range(flr_x(start_idx(j):stop_idx(j))));

            % higher values are low, so use minimum
            l_paw_max_y = cat(1, l_paw_max_y, min(fll_y(start_idx(j)):stop_idx(j)));
            r_paw_max_y = cat(1, r_paw_max_y, min(flr_y(start_idx(j)):stop_idx(j)));
            

            % interspersed with lick
            trial_length = length(lick_events);
            try
                lick_flag = any(lick_events(start_idx(j):stop_idx(j)+30));
            catch
                lick_flag = any(lick_events(start_idx(j):end));
            end

            if lick_flag
                lick_interspersed = cat(1, lick_interspersed, 1);
            else
                lick_interspersed = cat(1, lick_interspersed, 0);
            end

            % time between events

            % average speed
            fll_avg_speed = cat(1, fll_avg_speed, mean(fll_speed(start_idx(j):stop_idx(j))));
            flr_avg_speed = cat(1, flr_avg_speed, mean(flr_speed(start_idx(j):stop_idx(j))));

            % max speed
            fll_max_speed = cat(1, fll_max_speed, mean(fll_speed(start_idx(j):stop_idx(j))));
            flr_max_speed = cat(1, flr_max_speed, mean(flr_speed(start_idx(j):stop_idx(j))));

            % ratio of distance traveled
            fll_distance = cat(1, fll_distance, sum(fll_speed(start_idx(j):stop_idx(j))));
            flr_distance = cat(1, flr_distance, sum(flr_speed(start_idx(j):stop_idx(j))));

        end

        behavior_duration = cat(1,behavior_duration, (stop_idx - start_idx)./fs);
        behavior_label = cat(1, behavior_label, repmat(stroke_types(i), size(stop_idx,1), 1));
    end

end


y_diff = l_paw_max_y - r_paw_max_y;

%%

test = [behavior_duration, l_paw_xrange, l_paw_yrange, l_paw_start_y, ...
    l_paw_start_x, r_paw_xrange, r_paw_yrange, r_paw_start_y, ...
    r_paw_start_x, y_diff, lick_interspersed, fll_distance, flr_distance, ...
    imm_events, fll_avg_speed, flr_avg_speed, fll_max_speed, flr_max_speed];

[coeff, score, latent] = pca(test-mean(test));
[C, ia, ic] = unique(behavior_label);

%%
figure, hold on
% h=biplot(coeff(:,1:3),'scores',score(:,1:3));

for i = 1:size(behavior_label,1)
    switch behavior_label{i}
        case 'Elliptical'
            c = 'y';
            mkr = 'o';
        case 'Elliptical Asymmetric'
%             continue
            c = 'k';
            mkr = 'o';
        case 'Large Bilateral'
            continue
            c = 'g';
            mkr = 'o';
        case 'Left'
%             continue
            c = 'c';
            mkr = 'o';
        case 'Right'
%             continue
            c = 'm';
            mkr = 'o';
        case 'Left Asymmetric'
%             continue
            c = 'c';
            mkr = 'd';
        case 'Right Asymmetric'
%             continue
            c = 'm';
            mkr = 'd';
    end
    scatter3(score(i,1), score(i,2), score(i,3), c, mkr, 'filled', 'MarkerFaceAlpha', 0.25);
end


%%


figure, hold on
% h=biplot(coeff(:,1:3),'scores',score(:,1:3));

for i = 1:size(behavior_label,1)
    switch behavior_label{i}
        case 'Elliptical'
            c = 'k';
            mkr = 'o';
        case 'Elliptical Asymmetric'
            c = 'c';
            mkr = 'd';
        otherwise
            continue
    end
    scatter3(score(i,1), score(i,2), score(i,3), c, mkr)
end

%%
counter = 1;
for i = 1:length(h)
    if strcmpi(h(i).LineStyle, 'none')
        if counter <= ia(2)
            h(i).MarkerEdgeColor = 'k';
        elseif counter <= ia(3)
            h(i).MarkerEdgeColor = 'b';
        elseif counter <= ia(4)
            h(i).MarkerEdgeColor = 'm';
        elseif counter <= ia(5)
            h(i).MarkerEdgeColor = 'c';
        elseif counter <= ia(6)
            h(i).MarkerEdgeColor = 'g';
        else
            h(i). MarkerEdgeColor = 'r';
        end
        counter = counter + 1;

    end
end


%%











%%
events = read_boris('Y:\nick\behavior\grooming\1p\ECR2_thy1\20231116161015\ECR2_thy1_20231116161015_events.tsv');
figure, plot(events.("Right Asymmetric"))