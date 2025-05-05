clc, clear
addpath(genpath('Y:\nick\2p\code'))
load('atlas.mat')
load('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat');
[consolidated_file, experiment_path] = uigetfile('*.mat','Select consolidated neuron file.', 'Y:\nick\2p');
load([experiment_path, consolidated_file]);
% disp('[+] Finished loading consolidated neuron data')
% load([experiment_path, 'BN1_32puff_01_15fps_proc.mat']);
% disp('[+] Finished loading facemap data')


%%

imreg = imwarp(tforms{1}.linear_to_wfield.frame_wfi, tforms{1}.wfield_to_atlas.tform, 'interp', 'nearest','OutputView',imref2d(size(dorsalMaps.dorsalMapScaled)));




%%

x = [];
y = [];

for i = 1:length(F)
    F_smooth{i} = movmedian(F{i}, 3, 2);

    tmp = stat{i}(logical(iscell{i}(:,1)));
    for j = 1:length(tmp)
        x0 = double(tmp{j}.med(2));
        y0 = double(tmp{j}.med(1));        
        [x0, y0] = transformPointsForward(tforms{i}.roi_to_linear.tform, x0, y0);        
        [x0, y0] = transformPointsForward(tforms{i}.linear_to_wfield.tform, x0, y0);        
        [x0, y0] = transformPointsForward(tforms{i}.wfield_to_atlas.tform, x0, y0);
        x = [x; x0];
        y = [y; y0];
    end
    

end


figure, imshow(imreg, []), hold on
plot(x, y, 'r.', 'MarkerSize', 5)

for p = 1:length(dorsalMaps.edgeOutline)
    plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'w', 'LineWidth', 1.75);
end
set(gca, 'YDir', 'reverse');
title([num2str(length(x)), ' neurons'])

%%
[resampled_data, fs] = resample2maxFs(F_smooth, ops);
%% group left hemisphere
left_hem = 1:4;
right_hem = 5:6;
F_left = catcell(1, resampled_data(left_hem));
iscell_left = catcell(1, iscell(left_hem));
nloc_left = catcell(2, n_loc(left_hem));

good_idx = iscell_left(:,1) == 1;


[left_sorted, left_idx] = sort_neurons_by_region(F_left(good_idx,:), nloc_left(good_idx));

%% group right hemisphere

F_right = catcell(1, resampled_data(right_hem));
iscell_right = catcell(1, iscell(right_hem));
nloc_right = catcell(2,n_loc(right_hem));

good_idx = iscell_right(:,1) == 1;


[right_sorted, right_idx] = sort_neurons_by_region(F_right(good_idx,:), nloc_right(good_idx));

%%
sz = 10;
figure, subplot(sz, 1, 1:sz-3)
dataz = zscore([left_sorted; right_sorted], [], 2);
imagesc(dataz), axis off, caxis([-2 4])
colormap(flipud(colormap('gray')))
hold on
cols = {'c', 'm', 'y'};
ss = size(left_idx, 1);
left_tmp = left_idx;
for i = 1:3
    x = [1 size(left_sorted,2), size(left_sorted,2), 1];
    y = [find(left_idx == i, 1), ...
        find(left_idx == i, 1), ...
        find(left_idx == i+1, 1)-1, ...
        find(left_idx == i+1, 1)-1];
    patch(x,y,cols{i}, 'FaceAlpha', 0.12, 'EdgeColor', 'none')


    x = [1 size(left_sorted,2), size(left_sorted,2), 1];
    y = [ss + find(right_idx == i, 1), ...
        ss + find(right_idx == i, 1), ...
        ss + find(right_idx == i+1, 1)-1, ...
        ss + find(right_idx == i+1, 1)-1];
    patch(x,y,cols{i}, 'FaceAlpha', 0.12, 'EdgeColor', 'none')
end

%%

[U, s, V] = svd(dataz', 'econ');
V = s * V'; % multiply S into V, so only U and V from here on
U = U(:, 1:200); %reduce number of components
V = V(1:200, :); %reduce number of components


%%
t = xt(left_sorted, fs, 2);
sz = 10;
figure, subplot(sz, 1, 1:sz-3)
dataz = zscore(left_sorted, [], 2);
imagesc('Xdata', t, 'Ydata', 1:size(dataz,1), 'Cdata', dataz), axis off, caxis([-1 3]), axis tight
colormap(flipud(colormap('gray')))
hold on
cols = {'c', 'm', 'y'};
ss = size(left_idx, 1);
left_tmp = left_idx;
for i = 1:3
    x = [1 t(end), t(end), 1];
    y = [find(left_idx == i, 1), ...
        find(left_idx == i, 1), ...
        find(left_idx == i+1, 1)-1, ...
        find(left_idx == i+1, 1)-1];
    patch(x,y,cols{i}, 'FaceAlpha', 0.12, 'EdgeColor', 'none')

end


subplot(sz,1,sz-2)
plot(t, -U(:,1), 'k'), axis off, axis tight

t2 = xt(motion_1, 15);
subplot(sz,1,sz-1)
plot(t2, motion_1, 'k'), axis tight
ylim([min(motion_1), prctile(motion_1, 99.9)])
axis off

subplot(sz,1,sz)
plot(t2, motion_2, 'k'), axis tight
ylim([-prctile(motion_2, 99.9)/3, prctile(motion_2, 99.9)])
line([0.1 59.1], [-prctile(motion_2, 99.9)/5 -prctile(motion_2, 99.9)/5], 'Color', [0 0 0])
axis off

%%
% subplot(sz,1,sz)
% plot(pupil{1}.area(2:end), 'k'), axis tight
% axis off


