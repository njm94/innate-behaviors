

%%   do the same thing on ridge unique explained var

clear, clc

data_root = 'Y:\nick\behavior\grooming\1p';
mice = {'ECR2_thy1', 'GER2_ai94', 'HYL3_tTA', 'IBL2_tTA'};



for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    dff_path = [data_root, filesep, mice{j}, filesep, 'outputs'];
    h = openfig([dff_path, filesep, getAllFiles(dff_path, 'ridge_summary.fig')]);
    rightmove(:,:,j) = h.Children(2).Children.CData;
    leftmove(:,:,j) = h.Children(4).Children.CData;
    left(:,:,j) = h.Children(10).Children.CData;
    right(:,:,j) = h.Children(8).Children.CData;
    lick(:,:,j) = h.Children(6).Children.CData;
    elliptical(:,:,j) = h.Children(18).Children.CData;
    largeleft(:,:,j) = h.Children(16).Children.CData;
    largeright(:,:,j) = h.Children(14).Children.CData;
    bilateral(:,:,j) = h.Children(12).Children.CData;
    close(h)
end

%% overlay contours from diff mice
clc
nanmask = zeros(128, 128, length(mice));

for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    nanmask(:,:,j) = mask;
end
nanmask(nanmask==0) = nan;





%%
clc, clear a v
tmp = [];

thresh = 80;
figure, axis off, hold on
for j = 1:length(mice)

    vars = ["lick", "right", "left", "elliptical", "largeright", "largeleft", "bilateral"];
         
    for i = 1:length(vars)    
        test = eval(vars(i));
        test = test .* nanmask;
        test = test(:,:,j);

        v = prctile(test(:), thresh);
        a{i}(:,:,j) = test >= v;

        subplot(1,length(vars), i), axis off, hold on
            if v > 0
                contourf(flipud(test), [v v], 'FaceAlpha', 0.25)
    
                title(vars(i));
            end
%         end
    end
%     legend(labs, 'Location', 'Best')
end

%% compute pairwise dice
for i = 1:length(vars)
    count = 1;
    for j = 1:length(mice)
        for k = 1:length(mice)
            if k > j
                similarity(count, i) = dice(a{i}(:,:,j), a{i}(:,:,k));
                count = count + 1;
            end
        end
    end
end




[p,tbl,stats] = anova1(similarity);
[c,m,h,gnames] = multcompare(stats);

figure, boxplot(similarity, 'Colors', 'k'),
xticklabels(vars)
ylabel('Pairwise Dice Similarity Coefficient')
ax = gca;
ax.FontSize = 12;

title(['One-Way ANOVA, p = ', num2str(p, 2)])


%%




%%

selpath = uigetdir('Y:\nick\behavior\grooming\1p');

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


clc
visual = false;
cBetaRight = check_beta('Right', fullLabels, fullIdx, Umaster, fullBeta{1}, Vfull, mask, visual);
right = movmean(cBetaRight, 6, 3);

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