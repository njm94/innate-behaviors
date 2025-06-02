% Make behavior features for UMAP representation
%




%%

clc, clear
fileID = fopen('expt1_datalist.txt','r');
if isunix
    addpath('/home/user/Documents/grooming/utils')
else
    addpath('C:\Users\user\Documents\Nick\grooming\utils')
end

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);


mp_list = {'Y:\nick\behavior\grooming\2p\ETR2_thy1\20231113143925'; ...
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231113155903'; ...
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231115174148'; ...
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240729'; ...
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240731';
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240802';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240729';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240731';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240802';
    'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240729';
    'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240802'};
data_list{1} = [data_list{1}; mp_list];

current_mouse = '';

fs = 90 ;

thy1_idx = 1:7;
ai94_idx = 8:13;
camk_idx = 14:25;
mp_idx = 26:28;


%% 

% initialize feature and label matrices
bFilepath = {};
behavior_index = [];

behavior_duration = [];
behavior_label = [];

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


% Position features
fll_xpos_start = []; fll_ypos_start = []; % left paw (x,y) at start
flr_xpos_start = []; flr_ypos_start = []; % right paw (x,y) at start
fll_xpos_mid = []; fll_ypos_mid = []; % left paw (x,y) at midpoint
flr_xpos_mid = []; flr_ypos_mid = []; % right paw (x,y) at midpoint
fll_xpos_end = []; fll_ypos_end = []; % left paw (x,y) at end
flr_xpos_end = []; flr_ypos_end = []; % right paw (x,y) at end
fll_xpos_atmax_vL = []; fll_ypos_atmax_vL = []; % left paw (x,y) at time of max velocity of left paw
flr_xpos_atmax_vR = []; flr_ypos_atmax_vR = []; % right paw (x,y) at time of max velocity of right paw
fll_xpos_atmax_vR = []; fll_ypos_atmax_vR = []; % left paw (x,y) at time of max velocity of right paw
flr_xpos_atmax_vL = []; flr_ypos_atmax_vL = []; % right paw (x,y) at time of max velocity of left paw


fll_vx_atmax_vL = []; fll_vy_atmax_vL = [];
flr_vx_atmax_vR = []; flr_vy_atmax_vR = [];
fll_vx_atmax_vR = []; fll_vy_atmax_vR = [];
flr_vx_atmax_vL = []; flr_vy_atmax_vL = [];

corr_x = []; corr_y = [];

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
        dlc_vel_file = [data_dir, filesep, get_file_with_str(data_dir, 'vel.csv')];
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

    % remove lick, flail, and point events from event matrix
    idx = contains(events.Properties.VariableNames, 'Lick') | ...
        contains(events.Properties.VariableNames, 'Flail') | ...
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
    
    % get all position data relative to nose to account for variability
    % in camera positioning across trials
    flr_x = nose_x - dlc_pos(:,4);
    flr_y = nose_y - dlc_pos(:,5);    
    fll_x = nose_x - dlc_pos(:,7);
    fll_y = nose_y - dlc_pos(:,8);

    dlc_speed = readmatrix(dlc_speed_file);
    fll_speed = dlc_speed(:,1);
    flr_speed = dlc_speed(:,2);

    dlc_vel = readmatrix(dlc_vel_file);
    flr_vx = dlc_vel(:,4);
    flr_vy = dlc_vel(:,5);
    fll_vx = dlc_vel(:,7);
    fll_vy = dlc_vel(:,8);



    %% 

    stroke_types = stroke_events.Properties.VariableNames;

    for i = 1:length(stroke_types)
        idx = strcmpi(stroke_types{i}, boris.Behavior);
        behavior_frames = boris.ImageIndex(idx);
        
        % test = arr2idx(table2array(stroke_events(:,i)));
        % start_idx = test(:,1);
        % stop_idx = test(:,2);
        
        start_idx = behavior_frames(1:2:end);
        stop_idx = behavior_frames(2:2:end);

        behavior_index = cat(1, behavior_index, [start_idx, stop_idx]);

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


        fll_xpos_start = cat(1, fll_xpos_start, fll_x(start_idx));
        fll_ypos_start = cat(1, fll_ypos_start, fll_y(start_idx));
        flr_xpos_start = cat(1, flr_xpos_start, flr_x(start_idx));
        flr_ypos_start = cat(1, flr_ypos_start, flr_y(start_idx));

        fll_xpos_end = cat(1, fll_xpos_end, fll_x(stop_idx));
        fll_ypos_end = cat(1, fll_ypos_end, fll_y(stop_idx));
        flr_xpos_end = cat(1, flr_xpos_end, flr_x(stop_idx));
        flr_ypos_end = cat(1, flr_ypos_end, flr_y(stop_idx));

        for j = 1:length(start_idx)
            bFilepath = cat(1, bFilepath, boris_file);

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
                lick_flag = any(lick_events(start_idx(j)-3:stop_idx(j)+3));
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
            [maxfllspeed, argmaxfllspeed] = max(fll_speed(start_idx(j):stop_idx(j)));
            [maxflrspeed, argmaxflrspeed] = max(flr_speed(start_idx(j):stop_idx(j)));
            fll_max_speed = cat(1, fll_max_speed, maxfllspeed);
            flr_max_speed = cat(1, flr_max_speed, maxflrspeed);

            % position at max speed of same paw
            fll_xpos_atmax_vL = cat(1, fll_xpos_atmax_vL, fll_x(start_idx(j)+argmaxfllspeed));
            fll_ypos_atmax_vL = cat(1, fll_ypos_atmax_vL, fll_y(start_idx(j)+argmaxfllspeed));
            flr_xpos_atmax_vR = cat(1, flr_xpos_atmax_vR, flr_x(start_idx(j)+argmaxflrspeed));
            flr_ypos_atmax_vR = cat(1, flr_ypos_atmax_vR, flr_y(start_idx(j)+argmaxflrspeed));

            % position at max speed of opposite paw
            fll_xpos_atmax_vR = cat(1, fll_xpos_atmax_vR, fll_x(start_idx(j)+argmaxflrspeed));
            fll_ypos_atmax_vR = cat(1, fll_ypos_atmax_vR, fll_y(start_idx(j)+argmaxflrspeed));
            flr_xpos_atmax_vL = cat(1, flr_xpos_atmax_vL, flr_x(start_idx(j)+argmaxfllspeed));
            flr_ypos_atmax_vL = cat(1, flr_ypos_atmax_vL, flr_y(start_idx(j)+argmaxfllspeed));

            % midpoint position
            fll_xpos_mid = cat(1, fll_xpos_mid, fll_x(round(mean([start_idx(j), stop_idx(j)]))));
            fll_ypos_mid = cat(1, fll_ypos_mid, fll_y(round(mean([start_idx(j), stop_idx(j)]))));
            flr_xpos_mid = cat(1, flr_xpos_mid, flr_x(round(mean([start_idx(j), stop_idx(j)]))));
            flr_ypos_mid = cat(1, flr_ypos_mid, flr_y(round(mean([start_idx(j), stop_idx(j)]))));

            % velocity of paw at same paw's max speed (to capture direction)
            fll_vx_atmax_vL = cat(1, fll_vx_atmax_vL, fll_vx(start_idx(j)+argmaxfllspeed));
            fll_vy_atmax_vL = cat(1, fll_vy_atmax_vL, fll_vy(start_idx(j)+argmaxfllspeed));
            flr_vx_atmax_vR = cat(1, flr_vx_atmax_vR, flr_vx(start_idx(j)+argmaxflrspeed));
            flr_vy_atmax_vR = cat(1, flr_vy_atmax_vR, flr_vy(start_idx(j)+argmaxflrspeed));

            % velocity of paw at opposite paw's max speed (to capture direction)
            fll_vx_atmax_vR = cat(1, fll_vx_atmax_vR, fll_vx(start_idx(j)+argmaxflrspeed));
            fll_vy_atmax_vR = cat(1, fll_vy_atmax_vR, fll_vy(start_idx(j)+argmaxflrspeed));
            flr_vx_atmax_vL = cat(1, flr_vx_atmax_vL, flr_vx(start_idx(j)+argmaxfllspeed));
            flr_vy_atmax_vL = cat(1, flr_vy_atmax_vL, flr_vy(start_idx(j)+argmaxfllspeed));


            % distance traveled
            fll_distance = cat(1, fll_distance, sum(fll_speed(start_idx(j):stop_idx(j))));
            flr_distance = cat(1, flr_distance, sum(flr_speed(start_idx(j):stop_idx(j))));

            % correlation of paw trajectories
            corr_x = cat(1, corr_x, corr(fll_x, flr_x));
            corr_y = cat(1, corr_y, corr(fll_y, flr_y));

        end

        behavior_duration = cat(1,behavior_duration, (stop_idx - start_idx)./fs);
        behavior_label = cat(1, behavior_label, repmat(stroke_types(i), size(stop_idx,1), 1));
    end

end


%%

% use features that emphasize difference between right and left to
% better separate left/right actions

y_diff = l_paw_max_y - r_paw_max_y;
d_diff = fll_distance - flr_distance;
v_diff = fll_avg_speed - flr_avg_speed;
ystart_diff = fll_ypos_start - flr_ypos_start;
xstart_diff = fll_xpos_start - flr_xpos_start;
yrange_diff = l_paw_yrange - r_paw_yrange;
xrange_diff = l_paw_xrange - r_paw_xrange;
xmid_diff = fll_xpos_mid - flr_xpos_mid;
ymid_diff = fll_ypos_mid - flr_ypos_mid;
xend_diff = fll_xpos_end - flr_xpos_end;
yend_diff = fll_ypos_end - flr_ypos_end;
fll_range_ratio = l_paw_xrange./l_paw_yrange;
flr_range_ratio = r_paw_xrange./r_paw_yrange;

featureMat = zscore([behavior_duration, l_paw_xrange, l_paw_yrange, fll_ypos_start, ...
    fll_xpos_start, r_paw_xrange, r_paw_yrange, flr_ypos_start, ...
    flr_xpos_start, y_diff,  fll_distance, flr_distance, ...
    fll_avg_speed, flr_avg_speed, fll_max_speed, flr_max_speed, ...
    d_diff, v_diff, ystart_diff, xstart_diff, xrange_diff, yrange_diff, ...
    fll_xpos_atmax_vL, fll_ypos_atmax_vL, flr_xpos_atmax_vR, flr_ypos_atmax_vR, ...
    fll_xpos_mid, fll_ypos_mid, flr_xpos_mid, flr_ypos_mid, xmid_diff, ymid_diff, ...
    fll_xpos_end, fll_ypos_end, flr_xpos_end, flr_ypos_end, xend_diff, yend_diff, ...
    fll_xpos_atmax_vR, fll_ypos_atmax_vR, flr_xpos_atmax_vL, flr_ypos_atmax_vL, ... 
    fll_range_ratio, flr_range_ratio]);


%%
disp('Saving')
timenow = char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss'));
% save(['/media/user/teamshare/nick/behavior/grooming/', timenow, '_umap_test.mat'], 'featureMat', 'bFilepath', 'behavior_label', 'behavior_index')

