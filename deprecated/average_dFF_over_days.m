clc, clear


if ~isunix
    addpath('C:\Users\user\Documents\Nick\ridgeModel');
    addpath('C:\Users\user\Documents\Nick\ridgeModel\widefield')
    addpath('C:\Users\user\Documents\Nick\ridgeModel\smallStuff') 
    addpath('C:\Users\user\Documents\Nick\grooming\utils')
else
    addpath('/home/user/Documents/grooming/ridgeModel');
    addpath('/home/user/Documents/grooming/ridgeModel\widefield')
    addpath('/home/user/Documents/grooming/ridgeModel\smallStuff') 
    addpath('/home/user/Documents/grooming/utils')
end

fileID = fopen('expt1_datalist.txt','r');

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = '';

fs = 90 ;
[b, a] = butter(2, 0.01/(fs/2), 'high');

thy1_idx = 1:7;
ai94_idx = 8:13;
camk_idx = 14:25;
hyl3_idx = 14:19;
ibl2_idx = 20:25;

save_average_across_days = true;

load(fix_path('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat'));
%%
% spon_behavior_frames = cell(1,14);
% evoked_behavior_frames = cell(1,14);
bFrames = cell(1,14);
for j = thy1_idx %23:length(data_list{1})+1
    try
        data_dir = fix_path(data_list{1}{j});

        if isunix % working on linux computer - modify paths
            data_dir = strrep(data_dir, '\', '/');
            data_dir = strrep(data_dir, 'Y:/', '/media/user/teamshare/');
        end    
        disp(['Loading ' data_dir])
        % If there are multiple versions, use the most recent one
        b_frame_file = getAllFiles([data_dir, filesep, 'outputs'], 'behavior_frames.mat');
        if size(b_frame_file,1) > 1
            b_frame_file = sort(b_frame_file);
            b_frame_file = b_frame_file{end};
        end
        load([data_dir, filesep, 'outputs', filesep,b_frame_file])

        disp('Done loading')

        for i = 1:length(behavior_frames)
            if isempty(bFrames{i})
                bFrames{i} = behavior_frames{i};
            else
                bFrames{i} = cat(3, bFrames{i}, behavior_frames{i});
            end
        end

%         for i = 1:length(behavior_frames)
%             if j <= 2
%                 if isempty(spon_behavior_frames{i})
%                     spon_behavior_frames{i} = behavior_frames{i};
%                 else
%                     spon_behavior_frames{i} = cat(3, spon_behavior_frames{i}, behavior_frames{i});
%                 end
%             else
%                 if isempty(evoked_behavior_frames{i})
%                     evoked_behavior_frames{i} = behavior_frames{i};
%                 else, "bilateral"
%                     evoked_behavior_frames{i} = cat(3, evoked_behavior_frames{i}, behavior_frames{i});
%                 end
%             end
%         end

    catch
        new_mouse = true;
    end


end

%% 


figure
mouse_dir = fileparts(data_dir);
atlas_tform = load([mouse_dir, filesep, 'atlas_tform.mat']);
load([mouse_dir, filesep, 'mask.mat']);

for i = 1:length(bvars)

    subplot(4,4,i)
    if ~isempty(bFrames{i})
        mean_image = mask.*mean(bFrames{i},3)';
        mean_image = imwarp(mean_image, atlas_tform.tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
        imagesc(mean_image)
        title(bvars(i)), colorbar, xticks([]),  yticks([])
        hold on;
        for p = 1:length(dorsalMaps.edgeOutline)
            plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'w');
        end
        set(gca, 'YDir', 'reverse');
    end
    
end
savefig(gcf, [mouse_dir, filesep, 'outputs', filesep, char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '_dFF.fig'])


%%
figure
for i = 10:11%length(bvars)
    subplot(1,2,i-9)
    if ~isempty(evoked_behavior_frames{i})
        mean_image = mask.*mean(evoked_behavior_frames{i},3)';
    
        mean_image = imwarp(mean_image, atlas_tform.tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
        imagesc(mean_image)
        title(bvars(i)), colorbar, xticks([]),  yticks([])
        hold on;
        for p = 1:length(dorsalMaps.edgeOutline)
            plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'w');
        end
        set(gca, 'YDir', 'reverse');
    end
    
end


%%


figure
for i = 1:3%length(bvars)
    subplot(1,3,i)
    if ~isempty(behavior_frames{i})
        mean_image = mask.*mean(behavior_frames{i},3)';
    
        mean_image = imwarp(mean_image, atlas_tform.tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
        imagesc(mean_image)
        title(bvars(i)), colorbar, xticks([]),  yticks([])
        hold on;
        for p = 1:length(dorsalMaps.edgeOutline)
            plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'w');
        end
        set(gca, 'YDir', 'reverse');
    end
    
end


