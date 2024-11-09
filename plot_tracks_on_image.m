clear
addpath('C:\Users\user\Documents\Nick\ridgeModel');
addpath('C:\Users\user\Documents\Nick\ridgeModel\widefield')
addpath('C:\Users\user\Documents\Nick\ridgeModel\smallStuff') 
addpath('C:\Users\user\Documents\Nick\grooming\utils')



clc, clear
fileID = fopen('expt1_datalist.txt','r');

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = '';

fs = 90 ;

thy1_idx = 1:7;
ai94_idx = 8:13;
camk_idx = 14:25;
 

%%
for j = 19%:length(data_list{1})
    data_dir = data_list{1}{j};
    brain_file = [data_dir, filesep, getAllFiles(data_dir, 'cam0_svd')];
    beh_file = [data_dir, filesep, getAllFiles(data_dir, 'trim.mp4')];
    [mouse_root_dir, exp_date, ~] = fileparts(data_dir);
    [~, mouse_id, ~] = fileparts(mouse_root_dir);
    [expt_root_dir, ~, ~] = fileparts(mouse_root_dir);
    snippets_dir = [data_dir, filesep, 'snippets'];
    boris_file = [data_dir, filesep, getAllFiles(data_dir, '_events.tsv')];
    if ~isfile(boris_file), continue; end
    output_dir = [data_dir, filesep, 'outputs'];

%     if ~isempty(getAllFiles(output_dir, '.png'))
%         disp(['Tracks for ', data_dir, ' have been created already. Moving on...'])
%         continue
%     end

    if ~isfolder(output_dir)
        mkdir(output_dir)
    end

%     snippets_dir = [data_dir filesep 'snippets'];
    dlc_pos_file = [data_dir, filesep, getAllFiles(data_dir, '1030000.csv')];
    
%     [b, anno] = parse_snippets(snippets_dir);
    % read boris file, remove point event variables
    [events, b_idx, ~] = read_boris(boris_file);
    video_end = find(events.('Video End'));
    vid_end_idx = contains(events.Properties.VariableNames, 'Video End');
    events = removevars(events, vid_end_idx);
    b_idx(vid_end_idx) = [];
    % consolidate drop events
    drop_idx = contains(events.Properties.VariableNames, 'Drop');
    if any(drop_idx)
        events = removevars(events, drop_idx);
        b_idx(drop_idx) = [];
    end
    anno = events.Properties.VariableNames;




    % load DLC tracks
    disp('Loading DLC tracks')
    dlc_pos = readmatrix(dlc_pos_file);
    nose_x = median(dlc_pos(:,1));
    nose_y = median(dlc_pos(:,2));
    
    flr_x = dlc_pos(:,4);
    flr_y = dlc_pos(:,5);
    
    fll_x = dlc_pos(:,7);
    fll_y = dlc_pos(:,8);


    
    
    for jj = 1:length(b_idx)


        v = VideoReader(beh_file);
        img = read(v, round(mean(b_idx{jj}(1,:))));
    
        figure
%         imshow(imadjust(img, [], [], 0.5)), colormap gray, hold on, axis off, clim([500 5000])
        imshow(img(:,:,1), [0 200]), colormap gray, hold on, axis off, %clim([500 5000])
        
        for i = 1:size(b_idx{jj}, 1)
            plot(flr_x(b_idx{jj}(i,1):b_idx{jj}(i,2)), flr_y(b_idx{jj}(i,1):b_idx{jj}(i,2)), 'Color', [1 0 1 0.4])
            plot(fll_x(b_idx{jj}(i,1):b_idx{jj}(i,2)), fll_y(b_idx{jj}(i,1):b_idx{jj}(i,2)), 'Color', [0 1 1 0.4])
        end
        % center image on nose
        axis([nose_x - 150 nose_x + 150 nose_y-100 nose_y+200])
    
    
        fig = gcf;
        fig.Renderer = 'Painters';
        exportgraphics(gcf, [output_dir, filesep, char(anno{jj}),'.png'], 'Resolution', 300)
        close all
    
    end

end

%%
snippets_dir = [data_dir filesep 'snippets'];
dlc_pos_file = [data_dir, filesep, getAllFiles(data_dir, '1030000.csv')];

[b, anno] = parse_snippets(snippets_dir);

%%

% load DLC tracks
disp('Loading DLC tracks')
dlc_pos = readmatrix(dlc_pos_file);
flr_x = dlc_pos(:,4);
flr_y = dlc_pos(:,5);

fll_x = dlc_pos(:,7);
fll_y = dlc_pos(:,8);
%%

% figure, 
for j = 1:length(b)

    bvids = getAllFiles(snippets_dir, anno(j));
    if iscell(bvids)
        bvids = bvids{end};
    end
    img = loadtiff([snippets_dir filesep bvids]);
    img = img(:, :, 18);

    figure
    imshow(imadjust(img, [], [], 0.5)), colormap gray, hold on, axis off, clim([500 4000])
    
    for i = 1:min([size(b{j}, 1) 10])
        plot(flr_x(b{j}(i,1):b{j}(i,2)), flr_y(b{j}(i,1):b{j}(i,2)), 'Color', [1 0 1 0.75])
        plot(fll_x(b{j}(i,1):b{j}(i,2)), fll_y(b{j}(i,1):b{j}(i,2)), 'Color', [0 1 1 0.75])
        axis([200 460 20 270])
    end

fig = gcf;
fig.Renderer = 'Painters';
exportgraphics(gcf, ['Y:\nick\behavior\grooming\figures\', char(anno(j)),'.png'], 'Resolution', 300)

end

%%
