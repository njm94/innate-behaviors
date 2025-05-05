


clear, close all, clc
% cd('/media/user/teamshare/nick/behavior/grooming/code/')
% addpath('/media/user/teamshare/nick/behavior/grooming/code/ridgeModel')
% addpath('/media/user/teamshare/nick/behavior/grooming/code/ridgeModel/widefield')
% addpath('/media/user/teamshare/nick/behavior/grooming/code/ridgeModel/smallStuff')

addpath('C:\Users\user\Documents\Nick\ridgeModel');
addpath('C:\Users\user\Documents\Nick\ridgeModel\widefield')
addpath('C:\Users\user\Documents\Nick\ridgeModel\smallStuff') 

%%

fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = '';

fs = 90;
% [b, a] = butter(2, 0.01/(fs/2), 'high');
[b, a] = butter(1, [0.01 10]/(fs/2)); % bandpass

for j = 4%1:length(data_list{1})%1:length(data_list{1})+1
    % load experiment specific data into cell array
    ridge_dir = [data_list{1}{j}, filesep, 'ridge_outputs_video'];
    if ~isfolder(ridge_dir)
        disp('Ridge on Videp not performed yet for current experiment date')
        continue
    end
    
    
%     if any(contains(getAllFiles(ridge_dir), 'residuals'))
%         disp('Residual figure alreay created. skipping')
%         continue
%     end
    disp(['Starting ', data_list{1}{j}])

    disp('Loading ridge outputs')
    load([ridge_dir, filesep, getAllFiles(ridge_dir, 'cvFull')])
    data_dir = fileparts(data_list{1}{j});
    load([data_list{1}{j} filesep 'tform.mat'])
    load([data_dir filesep 'mask.mat'])
    snippets_dir = [data_list{1}{j}, filesep, 'snippets'];
    brain_file = [data_list{1}{j}, filesep, 'cam0_svd.mat'];
    dlc_speed_file = [data_list{1}{j}, filesep, getAllFiles(data_list{1}{j}, 'speed.csv')];

    
    
    disp('Loading brain data...')
    load(brain_file);
    % light-OFF signal may be captured in the brain data, resulting in a
    % massive filter artifact - remove a few frames just in case
    trial_length = size(V, 2) - 5; 
    Vbrain = s*V(:, 1:trial_length);
    
    Vbrain = filtfilt(b, a, Vbrain')';

    real_data = U*Vbrain;
%     Vfull = filtfilt(b,a,double(Vfull'))'; % filter 2nd time
    model_data = U*Vfull;
    
    disp('Calculating residuals')
    residuals = zscore(real_data, [], 2) - zscore(model_data,[],2);
    residuals = permute(reshape(residuals, 128, 128, []), [2 1 3]);
    residuals = evaluate_tform(residuals, tformEstimate);


    [behaviors, annotations, ~] = parse_snippets(snippets_dir);

    disp('Creating figure')

    nanmask = double(mask);
    nanmask(mask==0) = nan;

    disp('Loading DLC tracks')
    dlc_speed = readmatrix(dlc_speed_file);
    fll_speed = dlc_speed(1:trial_length,1);
    flr_speed = dlc_speed(1:trial_length,2);
    
    % binarize movement speed
    fll_speed = fll_speed > mean(fll_speed) + std(fll_speed);
    flr_speed = flr_speed >  mean(flr_speed) + std(flr_speed);
    
%     fll_speed = aggregate(fll_speed, 3);
%     flr_speed = aggregate(flr_speed, 3);
    fl_move = aggregate((fll_speed | flr_speed) , 3);

shgsh

    figure
    for i = 1:length(behaviors)
        test = [];
        if size(behaviors{1}, 1) > 0
            for jj = 1:size(behaviors{i},1)
                test = cat(3, test, residuals(:,:,behaviors{i}(jj,1):behaviors{i}(jj,2)));
            end
        end

        if ~isempty(test)
            subplot(3, 3, i)
            imagesc(nanmask.*mean(test,3))
            c=colorbar;
            c.Label.String = 'Residuals';
            title(annotations{i})
            xticks([])
            yticks([])
        end
    end
    
    disp('Saving figure')
    savefig(gcf, [ridge_dir, filesep, 'residuals',  char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')), '.fig'])

    clear U s V Vbrain real_data model_data residuals behaviors test

end



%%

testgroom = [];
for i = 1:size(behaviors{1},1)
    testgroom = cat(3, testgroom, residuals(:,:,behaviors{1}(i,1):behaviors{1}(i,2)));
end

testreach = [];
for i = 1:size(reach_idx,1)
    testreach = cat(3, testreach, residuals(:,:,reach_idx(i,1):reach_idx(i,2)));
end

testmove = [];
for i = 1:size(move_idx,1)
%     if move_idx(i,2)- move_idx(i,1) < 100
        testmove = cat(3, testmove, residuals(:,:,move_idx(i,1):move_idx(i,2)));
%     end
end


figure, imagesc([nanmask.*mean(testgroom,3), nanmask.*mean(testreach,3), nanmask.*mean(testmove,3)])
colorbar
caxis([-0.4 0.4])
colormap(bluewhitered)

%%

figure, imagesc(mean(testgroom,3))
[x,y] = ginput(1);

bdur = diff(move_idx,1,2);
[B, I] = sort(bdur);
newB = move_idx(I,:);
max_beh = max(bdur);
test = nan(size(move_idx,1), max_beh);
for i = 1:size(move_idx,1)
    test(i, 1:B(i)+1) = getTimeseries(residuals(:,:,newB(i,1):newB(i,2)), [x y], 2)';
end

figure,
imagesc(test), colorbar, caxis([-2 2]), 
colormap(bluewhitered)


%%


behavior_idx = zeros(1, trial_length);
for i = 1:length(behaviors)
    if isempty(behaviors{i}), continue; end
    for j = 1: size(behaviors{i},1)
        behavior_idx(behaviors{i}(j,1):behaviors{i}(j,2)) = 1;
    end
end
behavior_idx = aggregate(behavior_idx, 3);
behavior_idx = arr2idx(behavior_idx);

%%

groom_idx = zeros(1, trial_length);
for i = 1:4
    if isempty(behaviors{i}), continue; end
    for j = 1: size(behaviors{i},1)
        groom_idx(behaviors{i}(j,1):behaviors{i}(j,2)) = 1;
    end
end
groom_idx = aggregate(groom_idx, 3);
groom_idx = arr2idx(groom_idx);
%%

reach_idx = zeros(1, trial_length);
for i = 5:6
    if isempty(behaviors{i}), continue; end
    for j = 1: size(behaviors{i},1)
        reach_idx(behaviors{i}(j,1):behaviors{i}(j,2)) = 1;
    end
end
reach_idx = aggregate(reach_idx, 3);
reach_idx = arr2idx(reach_idx);

figure, plot(groom_idx)
hold on, plot(reach_idx)
%%
% mres = mean(residuals);
t = xt(mres, fs);
mres = squeeze(residuals(y,x,:));
figure, plot(t,mres, 'k')

cols = {'y', 'y', 'y', 'y', 'g', 'g', 'm'};
hold on
% for i = 1:length(behaviors)
%     patchplot(t(behaviors{i}), [min(mres) max(mres)], cols{i}, 0.25);
% end
i=7;
patchplot(t(behaviors{i}), [min(mres) max(mres)], cols{i}, 0.25);


%%



[s, f, t] = spectrogram(test,1000,100,1024, fs,'yaxis');

%%





%%


show_mov(residuals(:,:,100:200))


%%

rr = cell(1, length(behaviors));
for j = 1:length(behaviors)
    for i = 1:length(behaviors{j})
        rr{j} = [rr{j}, mres(behaviors{j}(i,1):behaviors{j}(i,2))];
    end
end


rrr = cellfun(@mean, rr)

%%
% load([selpath, filesep, 'outputs', filesep, 'cvFull.mat'])
load([selpath, filesep, 'mask.mat'])
% load([selpath, filesep, 'Umaster.mat'])
template = loadtiff([selpath, filesep, 'template.tif']);

disp('Done Loading')


%% average dff

% dff_fig = getAllFiles([selpath, filesep, 'outputs'], 'dFF_fig')
% h = openfig([selpath, filesep, 'outputs', filesep, '2024-03-08-07-56-29_dFF.fig']); % ecr2
% h = openfig([selpath, filesep, 'outputs', filesep, '2024-03-08-07-59-52_dFF.fig']); % ger2
% h = openfig([selpath, filesep, 'outputs', filesep, '2024-03-08-08-04-01_dFF.fig']); % hyl3
h = openfig([selpath, filesep, 'outputs', filesep, '2024-03-08-08-07-36_dFF.fig']); % ibl2



rightmove = h.Children(8).Children.CData;
leftmove = h.Children(10).Children.CData;
left = h.Children(16).Children.CData;
right = h.Children(14).Children.CData;
lick = h.Children(12).Children.CData;
elliptical = h.Children(24).Children.CData;
largeleft = h.Children(22).Children.CData;
largeright = h.Children(20).Children.CData;
dropright = h.Children(2).Children.CData;
close gcf



clc
thresh = 85;
labs = [];

vars = ["lick", "right", "elliptical", "rightmove"];
cols = {'m', 'y', 'c', 'k'};

figure, 
imagesc(double(template).*double(mask)), colormap gray, hold on, axis off

for i = 1:length(vars)

    test = eval(vars(i)).*mask;
    test(test==0) = nan;
    v = prctile(test(:), thresh);
    if v > 0
        contour(test, [v v], cols{i}, 'LineWidth', 2), hold on
        labs = [labs, vars(i)];
    end
end
legend(labs, 'Location', 'Best')


%% rige unique explained var
h = openfig([selpath, filesep, 'outputs', filesep, 'ridge_summary.fig']);
rightmove = h.Children(2).Children.CData;
leftmove = h.Children(4).Children.CData;
left = h.Children(10).Children.CData;
right = h.Children(8).Children.CData;
lick = h.Children(6).Children.CData;
elliptical = h.Children(18).Children.CData;
largeleft = h.Children(16).Children.CData;
largeright = h.Children(14).Children.CData;
dropright = h.Children(22).Children.CData;
close gcf

clc
thresh = 95;
labs = [];
% vars = ["leftmove", "rightmove", "left", "right", "lick", "largeright", "largeleft", "elliptical", "dropright"];
% cols = {'k', 'k', 'y', 'y', 'm', 'c', 'c', 'w', 'b'};

vars = ["left", "right"];
cols = {'c', 'm'};

figure, 
imagesc(double(template).*double(mask)), colormap gray, hold on, axis off

for i = 1:length(vars)

    test = eval(vars(i)).*mask;
    test(test==0) = nan;
    v = prctile(test(:), thresh);
    if v > 0
        contour(test, [v v], cols{i}, 'LineWidth', 2), hold on
        labs = [labs, vars(i)];
    end
end
legend(labs, 'Location', 'Best')

%%


% clc
visual = true;
cBetaRight = check_beta('Drop Hits Left', fullLabels, fullIdx, Ubrain, fullBeta{1}, Vfull, [], visual);
% right = movmean(cBetaRight, 6, 3);

%%
cBetaLeft = check_beta('Left', fullLabels, fullIdx, Umaster, fullBeta{1}, Vfull, mask, visual);
left = movmean(cBetaLeft, 6, 3);


cBetaElliptical = check_beta('Elliptical', fullLabels, fullIdx, Umaster, fullBeta{1}, Vfull, mask, visual);
elliptical = movmean(cBetaElliptical, 6, 3);


cBetaLL =  check_beta('LargeLeft', fullLabels, fullIdx, Umaster, fullBeta{1}, Vfull, mask, visual);
largeleft =  movmean(cBetaLL, 6, 3);

cBetaLR = check_beta('LargeRight', fullLabels, fullIdx, Umaster, fullBeta{1}, Vfull, mask, visual);
largeright = movmean(cBetaLR, 6, 3);



cBetaLick = check_beta('Lick', fullLabels, fullIdx, Umaster, fullBeta{1}, Vfull, mask, visual);
lick = movmean(cBetaLick, 6, 3);


cBetaAudio = check_beta('Audio', fullLabels, fullIdx, Umaster, fullBeta{1}, Vfull, mask, visual);
audio = movmean(cBetaAudio, 6, 3);

%%
fs = 90;
t0 = 45;
r = 22:180; % data range 45 is t0
t = xt(r, fs)-r(1)/(t0*2);

% figure, 
% subplot(1,2,1), plot(t,squeeze(movmean(cBetaRight(x(1), y(1), r), 6, 3))), hold on
% plot(t,squeeze(movmean(cBetaElliptical(x(1), y(1), r), 6, 3))),
% axis([t(1) t(end) -1 4])
% 
% subplot(1,2,2), plot(t,squeeze(movmean(cBetaRight(x(2), y(2), r), 6, 3))), hold on
% plot(t,squeeze(movmean(cBetaElliptical(x(2), y(2), r), 6, 3))),
% axis([t(1) t(end) -1 4])

%  plot(squeeze(largeright(x(1), y(1), 22:224)))

%%
r = 45:180;
thresh = 90;
labs = [];

figure, 
imagesc(template), colormap gray, hold on, axis off


test = mean(right(:,:,r), 3);
test(test==0) = nan;
v = prctile(test(:), thresh);
if v > 0
    contour(test, [v v], 'k', 'LineWidth', 2), hold on
    labs = [labs, "Right"];
end


test = mean(largeright(:,:,r), 3);
test(test==0) = nan;
v = prctile(test(:), thresh);
if v > 0
    contour(test, [v v], 'y', 'LineWidth', 2), hold on
    labs = [labs, "LargeRight"];
end


test = mean(elliptical(:,:,r), 3);
test(test==0) = nan;
v = prctile(test(:), thresh);
if v > 0
    contour(test, [v v], 'b', 'LineWidth', 2), hold on
    labs = [labs, "Elliptical"];
end

test = mean(lick(:,:,r), 3);
test(test==0) = nan;
v = prctile(test(:), thresh);
if v > 0
    contour(test, [v v], 'm', 'LineWidth', 2), hold on
    labs= [labs, "lick"];
end



legend(labs, 'Location', 'Best')

%%

thresh = 95;
test = mean(left, 3);
test(test==0) = nan;
v = prctile(test(:), thresh);
figure, contour(flipud(test), [v v], 'k', 'LineWidth', 2), hold on

test = mean(largeleft, 3);
test(test==0) = nan;
v = prctile(test(:), thresh);
contour(flipud(test), [v v], 'y', 'LineWidth', 2), hold on

test = mean(elliptical, 3);
test(test==0) = nan;
v = prctile(test(:), thresh);
contour(flipud(test), [v v], 'b', 'LineWidth', 2), hold on

test = mean(lick, 3);
test(test==0) = nan;
v = prctile(test(:), thresh);
contour(flipud(test), [v v], 'm', 'LineWidth', 2), hold on


legend({'Left', 'LargeLeft', 'Ellip', 'Lick'}, 'Location', 'Best')