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

%%

h = openfig('Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240729\arm1\roi_1\bin2x2x1\ROI_aligned.fig')
%%
ax2 = h.Children(2).Children.CData;
ax3 = h.Children(3).Children.CData;

%%

mask_ax2 = ax2>1;


se = strel('disk', 3);
boundary_ax2 = imdilate(edge(mask_ax2), se);



figure, imagesc(ax3), colormap gray
freezeColors()
green = cat(3, zeros(size(ax3)), ones(size(ax3)), zeros(size(ax3)));
hold on 
h = imshow(green); 
hold off 
set(h, 'AlphaData', 65535.*boundary_ax2) 

ax = gca;
exportgraphics(ax, fix_path(['Y:\nick\behavior\grooming\figures\', 'roi2linear_example.png']), 'Resolution', 300)

%%
close all
h = openfig('Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240729\arm1\linear_aligned.fig');

ax2 = h.Children(2).Children.CData;
ax3 = h.Children(3).Children.CData;


mask_ax2 = ax2>1;
se = strel('disk', 10);
mask_ax2 = imclose(mask_ax2, se);


se = strel('disk', 3);
boundary_ax2 = imdilate(edge(mask_ax2), se);



figure, imagesc(ax3), colormap gray
freezeColors()
green = cat(3, zeros(size(ax3)), ones(size(ax3)), zeros(size(ax3)));
hold on 
h = imshow(green); 
hold off 
set(h, 'AlphaData', 65535.*boundary_ax2) 

ax = gca;
exportgraphics(ax, fix_path(['Y:\nick\behavior\grooming\figures\', 'linear2widefield_example.png']), 'Resolution', 300)

