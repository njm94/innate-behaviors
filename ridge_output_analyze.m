clear, clc
selpath = uigetdir('Y:\nick\behavior\grooming\1p');

load([selpath, filesep, 'outputs', filesep, 'cvFull.mat'])
load([selpath, filesep, 'mask.mat'])
load([selpath, filesep, 'Umaster.mat'])
template = loadtiff([selpath, filesep, 'template.tif']);

disp('Done Loading')

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