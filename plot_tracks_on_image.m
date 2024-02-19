% clear
addpath C:\Users\user\Documents\Nick\grooming\utils
data_dir = 'Y:\nick\behavior\grooming\1p\HYL3_tTA\20231120043753';
data_dir = 'Y:\nick\behavior\grooming\1p\HYL3_tTA\20231116184954';
% data_dir = 'Y:\nick\behavior\grooming\1p\HYL3_tTA\20231112090820';
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
