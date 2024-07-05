% 

clear, clc
%
addpath('C:\Users\user\Documents\Nick\grooming\utils')
data_path = 'Y:\pankaj\water_reaching_project\data\matlab_outputs';

reaching_files = {'HR1_2019_4_12data_compile.mat', 'HR2_2019_4_12data_compile.mat', ...
    'IJ2_2020_7_27data_compile.mat', 'IJ3_2020_7_27data_compile.mat', ...
    'IK1_2019_4_12data_compile.mat', 'IK2_2019_4_12data_compile.mat', ...
    'IO1_2020_7_27data_compile.mat', 'IQ1_2020_7_27data_compile.mat', ...
    'IQ2_2020_7_27data_compile.mat'};

% reaching_files = {'HR1_2019_4_12data_compile.mat', 'HR2_2019_4_12data_compile.mat', ...
%     'IK1_2019_4_12data_compile.mat', 'IK2_2019_4_12data_compile.mat'};

fs=60;
%%

for j = 6%length(reaching_files)
    disp(['Running ', data_path, filesep, reaching_files{j}])
    load([data_path, filesep, reaching_files{j}], 'dFF', 'rewardTime', ...
        'trialResult', 'lickTime', 'paw', 'beh', 'whisker_time', 'fluo_frame', ...
        'dff_fluo')
    N = length(rewardTime);

    figure, imagesc(beh{1}(:,:,1)), hold on, colormap gray
    for i = 1:length(trialResult)/2
        if strcmpi(trialResult{i}, 'Success')
             plot(paw{i}(:,1), paw{i}(:,2)), hold on, plot(paw{i}(:,3), paw{i}(:,4))
        end
    end
    
    r = drawpolygon(LineWidth=3,Color="y");
    l = drawpolygon(LineWidth=3,Color='m');

    mask = draw_roi(fluo_frame, 2);
    mask = double(mask);
    mask(mask==0) = nan;


    figure, hold on
    reach_response_nohemo = [];
    lick_response_nohemo = [];
    reach_response_hemo_corrected = [];
    lick_response_hemo_corrected = [];
    reach_response_non_overlapping = [];
    lick_response_non_overlapping = [];
    whisker_response = [];

    for i = 1:N
        if strcmpi(trialResult{i}, 'Success')

            lick_idx = inROI(l, paw{i}(:,1), paw{i}(:,2));
            reach_idx = inROI(r, paw{i}(:,1), paw{i}(:,2));

            % remove any events before reward delivered
            reach_idx(1:round(rewardTime(i)*fs)) = 0;
            lick_idx(1:round(rewardTime(i)*fs)) = 0;

            % smooth lick and reach indices with morphological opening
            w = round(fs/3);
            lick_idx = movmax(movmin(lick_idx, w), w);
            reach_idx = movmax(movmin(reach_idx, w), w);
            
            % exclude overlapping period from both signals
            overlap = lick_idx & reach_idx;
            reach_idx_non_overlapping = reach_idx;
            reach_idx_non_overlapping(overlap) = 0;

            lick_idx_non_overlapping = lick_idx;
            lick_idx_non_overlapping(overlap) = 0;
            
            plot(xt(reach_idx,fs), reach_idx-i, 'k')
            plot(xt(lick_idx,fs), lick_idx-i, 'm')
            plot(rewardTime(i), -i, 'o')
            plot(xt(whisker_time{i}, fs), whisker_time{i}-i, 'k')
            
            reach_response_nohemo = cat(3, reach_response_nohemo, dff_fluo{i}(:,:,reach_idx(1:min([length(reach_idx), size(dff_fluo{i},3)]))));
            lick_response_nohemo = cat(3, lick_response_nohemo, dff_fluo{i}(:,:,lick_idx(1:min([length(lick_idx), size(dff_fluo{i},3)]))));

            reach_response_hemo_corrected = cat(3, reach_response_hemo_corrected, dFF{i}(:,:,reach_idx(1:min([length(reach_idx), size(dFF{i},3)]))));
            lick_response_hemo_corrected = cat(3, lick_response_hemo_corrected, dFF{i}(:,:,lick_idx(1:min([length(lick_idx), size(dFF{i},3)]))));
            reach_response_non_overlapping = cat(3, reach_response_non_overlapping, dFF{i}(:,:,reach_idx_non_overlapping(1:min([length(reach_idx_non_overlapping), size(dFF{i},3)]))));
            lick_response_non_overlapping = cat(3, lick_response_non_overlapping, dFF{i}(:,:,lick_idx_non_overlapping(1:min([length(lick_idx_non_overlapping), size(dFF{i},3)]))));
            whisker_response = cat(3, whisker_response, dFF{i}(:,:,whisker_time{i}));            
        end
    end

    disp('Saving')
    save([data_path, filesep, reaching_files{j}(1:end-4), '_processed_', char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '.mat'], ...
        'reach*', 'lick*', 'r', 'l', 'mask', 'whisker*')
end

%%
clear, clc
addpath('C:\Users\user\Documents\Nick\grooming\utils')
data_path = 'Y:\pankaj\water_reaching_project\data\matlab_outputs';

% reaching_files_to_analyze = {'HR1_2019_4_12data_compile_processed_2024-06-20-11-37-34.mat',...
%     'HR2_2019_4_12data_compile_processed_2024-06-20-11-51-05.mat', ...
%     'IK1_2019_4_12data_compile_processed_2024-06-20-11-58-33.mat', ...
%     'IK2_2019_4_12data_compile_processed_2024-06-20-12-05-24.mat'};

reaching_files_to_analyze = {'HR1_2019_4_12data_compile_processed_2024-06-27-17-14-32.mat',...
    'HR2_2019_4_12data_compile_processed_2024-06-27-17-24-25.mat', ...
    'IK1_2019_4_12data_compile_processed_2024-06-27-17-30-00.mat', ...
    'IK2_2019_4_12data_compile_processed_2024-06-27-17-37-26.mat'};

thresh = 85;
fs = 60;

load('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat');


figure
for i = 1:length(reaching_files_to_analyze)
    disp(['Loading ', reaching_files_to_analyze{i}])
    load([data_path, filesep, reaching_files_to_analyze{i}])


    % avg_reaching_map = mean(test(:,:,prewin+1:prewin+1+round(0.5*fs),:), [4 3]).*mask;
    avg_reaching_map = mean(reach_response_nohemo(:,:,1:round(size(reach_response_nohemo,3)/3)),3) .* mask;
    avg_reaching_map = imwarp(avg_reaching_map, tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
    avg_reaching_map(avg_reaching_map==0) = nan;
    
    if i == 2
        subplot(1,2,1)
        imagesc(avg_reaching_map); colorbar, hold on
        for p = 1:length(dorsalMaps.edgeOutline)
            plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'w');
        end
        xticks([])
        yticks([])
    end

    subplot(1,2,2)
    axis off, hold on
    vr = prctile(avg_reaching_map(:), thresh);
    contourf(avg_reaching_map, [vr vr], 'FaceAlpha', 0.25)
    for p = 1:length(dorsalMaps.edgeOutline)
        plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k');
    end
    set(gca, 'YDir', 'reverse');

    disp('Halfway done')
end


%%


figure
for i = 1:length(reaching_files_to_analyze)
    disp(['Loading ', reaching_files_to_analyze{i}])
    load([data_path, filesep, reaching_files_to_analyze{i}])


    % avg_reaching_map = mean(test(:,:,prewin+1:prewin+1+round(0.5*fs),:), [4 3]).*mask;
    avg_reaching_map = mean(reach_response_hemo_corrected,3) .* mask;
%     avg_reaching_map = imwarp(avg_reaching_map, tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
%     avg_reaching_map(avg_reaching_map==0) = nan;
    

        subplot(1,4,i)
        imagesc(avg_reaching_map); colorbar, hold on
%         for p = 1:length(dorsalMaps.edgeOutline)
%             plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'w');
%         end
        xticks([])
        yticks([])

end