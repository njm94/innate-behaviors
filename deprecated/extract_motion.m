% Extract bouts of non-grooming motion

clear, close all, clc
%%
fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
for j = 4%1:length(data_list{1})
    data_dir = data_list{1}{j};
    disp(['Starting ' data_dir])
    fPath = [data_dir filesep 'ridge_outputs' filesep];
    if ~isdir(fPath), mkdir(fPath); end
%     if isfile([fPath, 'summary.fig'])
%         disp([data_dir, ' already processed. Skipping...'])
%         continue
%     end

    
    exp_date = strfind(data_dir, filesep);
    exp_date = data_dir(exp_date(end)+1:end);
    fs = 90;
    
    beh_file = [data_dir, filesep, get_file_with_str(data_dir, [exp_date, '_svd'])];
    ME_file = [data_dir, filesep, get_file_with_str(data_dir, 'MEsvd')];
    frame_file = [data_dir, filesep, get_file_with_str(data_dir, 'singleFrame')];
    grooming_file = [data_dir, filesep, 'grooming_events_roi_filtered.mat'];
    dlc_pos_file = [data_dir, filesep, get_file_with_str(data_dir, '1030000.csv')];
    dlc_speed_file = [data_dir, filesep, get_file_with_str(data_dir, 'speed.csv')];

    if any([~isfile(dlc_pos_file), ~isfile(beh_file), ~isfile(ME_file), ~isfile(grooming_file)])
        disp('Missing one or more of the required files. Skipping...')
        continue
    end

    disp('Loading DLC tracks')
    dlc_speed = readmatrix(dlc_speed_file);
    fll_speed = dlc_speed(:,1);
    flr_speed = dlc_speed(:,2);
    dlc_pos = readmatrix(dlc_pos_file);
    fll_x = dlc_pos(:,7);
    fll_y = dlc_pos(:,8);
    flr_x = dlc_pos(:,4);
    flr_y = dlc_pos(:,5);

    %% load grooming events
    disp('Loading grooming events')
    load(grooming_file)
    
    % separate ipsilateral, contralateral, and bilateral events
    ipsi = zeros(size(FL_L));
    contra = zeros(size(FL_L));
    bilat = zeros(size(FL_L));
    for i = 1:size(event_type, 1)
        if contains(event_type(i,:), 'ipsi')
            ipsi(all_events(i,1):all_events(i,2)) = 1;
        elseif contains(event_type(i,:), 'cont')
            contra(all_events(i,1):all_events(i,2)) = 1;
        elseif contains(event_type(i,:), 'bi')
            bilat(all_events(i,1):all_events(i,2)) = 1;
        else
            disp('wrong string in event type')
        end
    end

    %% try to find individual strokes using autocorrelation
    figure, 
    for i = 1:size(all_events,1)
        subplot(2,1,1)
        left_x = fll_x(all_events(i,1):all_events(i,2));
        left_y = fll_y(all_events(i,1):all_events(i,2));
        right_x = flr_x(all_events(i,1):all_events(i,2));
        right_y = flr_y(all_events(i,1):all_events(i,2));
        plot(right_x, -right_y)
        hold on
        plot(left_x, -left_y)
        title(event_type(i,:))
        
        hold off
        subplot(2,1,2)
%         [acf_x, lags_x] = xcorr(right_x);
        acf_x = cconv(right_x,conj(fliplr(right_x)),length(right_x));
        acf_y = cconv(right_y,conj(fliplr(right_y)),length(right_y));
%         [acf_y, lags_y] = xcorr(right_y);
        plot(zscore(acf_x))
        hold on
        plot(zscore(acf_y))
        hold off
        h = waitforbuttonpress;

    end


    %%
    thresh_left = mean(fll_speed) + std(fll_speed);
    thresh_right = mean(flr_speed) + std(flr_speed);
    fll_movement = fll_speed > thresh_left;
    flr_movement = flr_speed > thresh_right;

    fll_movement = movmax(fll_movement, 35);
    fll_movement = fll_movement & ~(ipsi | contra | bilat)';
    flr_movement = movmax(flr_movement, 35);
    flr_movement = flr_movement & ~(ipsi | contra | bilat)';
    


%     fll_movement


%%

figure,
t = xt(fll_speed, fs);
plot(t, fll_speed), hold on
plot(t, flr_speed)
axis tight; yL = ylim;
patchplot()
% patchplot(t(all_events(contains(cellstr(event_type), 'ipsi'), :)), yL, 'c', 0.3)
% patchplot(t(all_events(contains(cellstr(event_type), 'contra'), :)), yL, 'm', 0.3)
% patchplot(t(all_events(contains(cellstr(event_type), 'bi'), :)), yL, 'y', 0.3)

%%
   
end
  %%

function f = get_file_with_str(data_dir, str_in)
file_list = getAllFiles(data_dir);
file_idx = contains(file_list, str_in);
f = file_list{file_idx};
end
