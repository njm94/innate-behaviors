% 

clear, clc
%%
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

for j = 4%length(reaching_files)
    disp(['Running ', data_path, filesep, reaching_files{j}])
    load([data_path, filesep, reaching_files{j}], 'dFF', 'rewardTime', ...
        'trialResult', 'lickTime', 'paw', 'beh', 'whisker_time', 'fluo_frame')
    N = length(rewardTime);

    figure, imagesc(beh{1}(:,:,1)), hold on, colormap gray
    for i = 1:length(trialResult)/2
        if strcmpi(trialResult{i}, 'Success')
             plot(paw{i}(:,1), paw{i}(:,2)), hold on, plot(paw{i}(:,3), paw{i}(:,4))
        end
    end
    
    r = drawpolygon(LineWidth=3,Color="y");
    l = drawpolygon(LineWidth=3,Color='m');

    mask = draw_roi(fluo_frame, 4);
    mask = double(mask);
    mask(mask==0) = nan;


    figure, hold on
    reach_response = [];
    lick_response = [];
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
            reach_idx(overlap) = 0;
            lick_idx(overlap) = 0;
            
            plot(xt(reach_idx,fs), reach_idx-i, 'k')
            plot(xt(lick_idx,fs), lick_idx-i, 'm')
            plot(rewardTime(i), -i, 'o')
            plot(xt(whisker_time{i}, fs), whisker_time{i}-i, 'k')

            reach_response = cat(3, reach_response, dFF{i}(:,:,reach_idx(1:min([length(reach_idx), size(dFF{i},3)]))));
            lick_response = cat(3, lick_response, dFF{i}(:,:,lick_idx(1:min([length(lick_idx), size(dFF{i},3)]))));
            whisker_response = cat(3, whisker_response, dFF{i}(:,:,whisker_time{i}));            
        end
    end

    save([data_path, filesep, reaching_files{j}(1:end-4), '_processed_', char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '.mat'], ...
        'reach*', 'lick*', 'r', 'l', 'mask', 'whisker*')
end

%%



%%
clc
for j = 1%:length(reaching_files)
    disp(['Running ', data_path, filesep, reaching_files{j}])
    load([data_path, filesep, reaching_files{j}], 'dFF', 'rewardTime', ...
        'trialResult', 'lickTime', 'paw', 'beh', 'fluo_frame')
    N = length(rewardTime);
    clear reach_response lick_response

    figure, imshow(beh{1}(:,:,:,100)), hold on
    roi_flag = true;
    count = 1;
    for i = 1:N
        if (rewardTime(i) > 0) && strcmpi(trialResult{i}, 'Success')
            tlen = size(dFF{i},3);

            plot(paw{i}(:,1), paw{i}(:,2)), hold on, plot(paw{i}(:,3), paw{i}(:,4))
            
            x = paw{i}(1:tlen,1);
            y = paw{i}(1:tlen,2);

            for k = 1:length(x)
                lick_idx(k) = lick_roi(y(k), x(k), 1);
                reach_idx(k) = reach_roi(y(k), x(k), 1);
            end

            reach_response(:,:,count) = mean(dFF{i}(:,:,reach_idx), 3);
            lick_response(:,:,count) = mean(dFF{i}(:, :, lick_idx), 3);
            count = count +1;
            

%             disp(round(lickTime(i)*fs) - round(rewardTime(i)*fs))
%             try
%                 reach_response(:,:,count) = mean(dFF{i}(:,:,round(rewardTime(i)*fs):round(lickTime(i)*fs)),3);
%                 lick_response(:,:,count) = mean(dFF{i}(:,:,round(lickTime(i)*fs):end), 3);
%                 count = count + 1;
%             catch
%                 disp('')
%                 continue
%             end
        end
    end
    avg_reaching_map(:,:,j) = mean(reach_response, 3);
    avg_licking_map(:,:,j) = mean(lick_response, 3);
end


%%



 count = 1;
    for i = 1:N
        if (rewardTime(i) > 0) && strcmpi(trialResult{i}, 'Success')
            tlen = size(dFF{i},3);

            plot(paw{i}(:,1), paw{i}(:,2)), hold on, plot(paw{i}(:,3), paw{i}(:,4))
            
            x = round(paw{i}(1:tlen,1));
            y = round(paw{i}(1:tlen,2));

            clear lick_idx reach_idx

            for k = 1:length(x)
                lick_idx(k) = lick_roi(y(k), x(k), 1);
                reach_idx(k) = reach_roi(y(k), x(k), 1);
            end

            reach_response(:,:,count) = mean(dFF{i}(:,:, logical(reach_idx)), 3);
            lick_response(:,:,count) = mean(dFF{i}(:, :, logical(lick_idx)), 3);
            count = count +1;
            
        end
    end



%%


mask = draw_roi(fluo_frame, 2);
mask(mask==0) = nan;
%%

thresh = 80;
% avg_reaching_map = mean(test(:,:,prewin+1:prewin+1+round(0.5*fs),:), [4 3]).*mask;
avg_reaching_map = mean(reach_response,3) .* mask;
avg_licking_map =mean(lick_response,3) .* mask;
vr = prctile(avg_reaching_map(:), thresh);
vl = prctile(avg_licking_map(:), thresh);
figure, axis off, hold on
contourf(flipud(avg_reaching_map), [vr vr], 'FaceAlpha', 0.25)
contourf(flipud(avg_licking_map), [vl vl], 'FaceAlpha', 0.25)

