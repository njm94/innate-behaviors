% Frequency analysis of grooming behaviors to identify automatic vs
% voluntary behaviors
%
% Expected outcome - Automatic behaviors occur with stereotyped frequencies
%



%% average events

clc, clear


if ~isunix
    addpath('C:\Users\user\Documents\Nick\grooming\utils')
else
    addpath('/home/user/Documents/grooming/utils')
end




fileID = fopen('expt1_datalist.txt','r');

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = '';

fs = 90 ; % sampling rate


thy1_idx = 1:7;
ai94_idx = 8:13;
camk_idx = 14:25;

save_average_across_days = true;

%%
 
for j = camk_idx%1:length(data_list{1})+1

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
%%
            if save_average_across_days
            end
            
            %%
            % all mice completed - break the loop
            if j == length(data_list{1})+1, break; end
        end        
        count = 1;
        current_mouse = mouse_id;

        fPath = [mouse_root_dir filesep 'outputs' filesep];

    end

    try
        boris_file = [data_dir, filesep, get_file_with_str(data_dir, 'events.tsv')];
        dlc_pos_file = [data_dir, filesep, get_file_with_str(data_dir, '1030000.csv')];
        dlc_speed_file = [data_dir, filesep, get_file_with_str(data_dir, 'speed.csv')];
    catch
        disp('Missing one or more of the required files. Skipping...')
        continue
    end

    %% load grooming events from BORIS file
    events = read_boris(boris_file);

    % create the aggregated behavior episodes - if none exceed 15s, skip
    % consolidate lick events
    lick_idx = contains(events.Properties.VariableNames, 'Lick');
    lick_events = events(:,lick_idx);
    lick_events = any(table2array(lick_events),2);
    
    % remove lick and point events from event matrix
    rmidx = contains(events.Properties.VariableNames, 'Lick') | ...
        contains(events.Properties.VariableNames, 'Drop') | ...
        contains(events.Properties.VariableNames, 'Video') | ...
        contains(events.Properties.VariableNames, 'Flail') ;
    stroke_events = removevars(events, rmidx);
    labels = stroke_events.Properties.VariableNames;
    clear stroke_on
    for ii = 1:size(stroke_events,2)
        stroke_on(:,ii) = get_on_time(table2array(stroke_events(:,ii)));
    end

    bmat = single(any(stroke_on,2));
    brate = movsum(bmat, fs/2)*2;

    % define all automatic grooming episodes as those which exist within an
    % epsiode that at one point has a rate which exceeds THRESH behaviors 
    % per second.

    thresh = 3;

    [~, locs] = findpeaks(brate, 'MinPeakHeight', thresh+0.1, 'MinPeakProminence', thresh+1);
    % get start and end indices for each episode
    for ii = 1:size(locs,1)
        idx(ii,1) = locs(ii) - find(~flipud(brate(1:locs(ii))), 1);
        idx(ii,2) = locs(ii) + find(~brate(locs(ii):end), 1);
    end
    idx = unique(idx, 'rows'); % remove repeats
    
    % load DLC tracks
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


    % auto_count = zeros(1,size(stroke_on,2));
    % stroke_on_copy = stroke_on;
    % for ii = 1:size(idx,1)
    %     auto_count = auto_count + sum(stroke_on(idx(ii,1):idx(ii,2),:));
    %     stroke_on_copy(idx(ii,1):idx(ii,2), :) = 0;
    % end
    % vol_count = sum(stroke_on_copy);


    figure
    auto_idx = idx2arr(idx, size(stroke_events,1))';
    for i = 1:size(stroke_events,2)
        tmp = logical(table2array(stroke_events(:,i))) & auto_idx;
        
        subplot(2,7,i), hold on
        plot(flr_x(tmp), -flr_y(tmp))
        plot(fll_x(tmp), -fll_y(tmp))
        axis([200 450 -250 0])
        title(stroke_events.Properties.VariableNames{i})
        
        tmp = logical(table2array(stroke_events(:,i))) & ~auto_idx;
        subplot(2,7,i+7), hold on
        plot(flr_x(tmp), -flr_y(tmp))
        plot(fll_x(tmp), -fll_y(tmp))
        axis([200 450 -250 0])
        % title(stroke_events.Properties.VariableNames{i})
    end


end

