clear
addpath C:\Users\user\Documents\Nick\grooming\utils
data_dir = 'Y:\nick\behavior\grooming\1p\HYL3_tTA\20231120043753';
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

figure, 
for j = 1:length(b)

    bvids = getAllFiles(snippets_dir, anno(j));
    img = loadtiff([snippets_dir filesep bvids{end}]);
    img = img(:, :, 16);

    subplot(3,3,j)
    imagesc(imadjust(img)), colormap gray, hold on, axis off, clim([1000 40000])
    
    for i = 1:min([size(b{j}, 1) 10])
        plot(flr_x(b{j}(i,1):b{j}(i,2)), flr_y(b{j}(i,1):b{j}(i,2)), 'Color', [1 0 1 0.5])
        plot(fll_x(b{j}(i,1):b{j}(i,2)), fll_y(b{j}(i,1):b{j}(i,2)), 'Color', [0 1 1 0.5])
        axis([200 460 20 250])
    end
    title(anno(j))
end