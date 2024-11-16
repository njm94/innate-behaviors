% plot forelimb and hindlimb responses
clear
addpath(fix_path('C:\Users\user\Documents\Nick\grooming\utils'))
hl_response = readmatrix(fix_path('Y:\nick\behavior\grooming\2p\IDR3_tTA6s\dalsa\hl_response.csv'));
fl_response = readmatrix(fix_path('Y:\nick\behavior\grooming\2p\IDR3_tTA6s\dalsa\fl_response.csv'));
fs = 40;
t = xt(hl_response, fs, 1);

close all
figure, 
plot(t-3,hl_response(:,2), 'k', 'LineWidth', 2)
axis([-3 3 -5 20])
hold on
vline(0, 'k--')
xlabel('Time (s)')
ylabel('\DeltaF/F_0 (%)')
ax = gca;
saveas(ax, fix_path(['Y:\nick\behavior\grooming\figures\', 'HL_response', '.svg']))


figure, 
plot(t-3,fl_response(:,2), 'k', 'LineWidth', 2)
hold on
axis([-3 3 -2 10])
vline(0, 'k--')
yrange = ylim;
xlabel('Time (s)')
ylabel('\DeltaF/F_0 (%)')

ax = gca;
saveas(ax, fix_path(['Y:\nick\behavior\grooming\figures\', 'FL_response', '.svg']))

%%


load(fix_path('C:\Users\user\Documents\Nick\grooming\utils\allenDorsalMap.mat'))
figure
axis equal off; hold on;
for p = 1:length(dorsalMaps.edgeOutline)
    plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k', 'LineWidth', 2);
end
set(gca, 'YDir', 'reverse');
ax = gca;
saveas(ax, fix_path(['Y:\nick\behavior\grooming\figures\', 'AllenAtlas', '.svg']))
