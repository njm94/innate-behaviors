clear, clc
[wfi_file, wfi_path] = uigetfile('*.tif','Select widefield reference image.', 'Y:\nick\behavior\grooming\2p');
data_wfi = loadtiff([wfi_path, wfi_file]);
frame_wfi = imadjust(uint16(imflatfield(mean(data_wfi, 3), 100)));

[comp_file, comp_path] = uigetfile('*.tif','Select sensory composite image.', 'Y:\nick\behavior\grooming\2p');
data_comp = loadtiff([comp_path, comp_file]);
frame_comp = imadjust(uint16(data_comp(:,:,1)));
%%

[optimizer,metric] = imregconfig("monomodal");
tform = imregtform(frame_wfi, frame_comp, 'rigid', optimizer, metric);

movingReg = imwarp(frame_wfi, tform,"OutputView",imref2d(size(frame_comp)));

%%

save([wfi_path, 'template2composite_tform.mat'], 'tform', "frame_wfi", "frame_lin", "moving_registered");
close all;


%%
load('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat');

%% Do registration and transformation on atlas

tform = align_recording_to_allen(imadjust(frame_wfi), {}, 1); % align <-- input any function of data here
invT=pinv(tform.T); % invert the transformation matrix
invT(1,3)=0; invT(2,3)=0; invT(3,3)=1; % set 3rd dimension of rotation artificially to 0
invtform=tform; invtform.T=invT; % create the transformation with invtform as the transformation matrix
load('atlas.mat')
% atlas=imwarp(atlas,invtform, 'interp','nearest', 'OutputView', imref2d(size(frame_wfi))); % setting the 'OutputView' is important
% atlas=round(atlas);

%% Save outputs

imagereg = imwarp(imadjust(frame_wfi),tform,'interp','nearest','OutputView',imref2d(size(dorsalMaps.dorsalMapScaled)));

% widefield image with atlas overlaid
figure
% imagesc(imflatfield(imagereg, 100)); colormap gray
imagesc(imagereg); colormap gray
axis equal off; hold on;
for p = 1:length(dorsalMaps.edgeOutline)
    plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1));
end
set(gca, 'YDir', 'reverse');
saveas(gcf, [wfi_path, 'atlas_aligned.fig'])


%%

% Warped atlas image
figure; subplot(1,2,1); imagesc(imflatfield(frame_wfi,100)); axis off, grid on, colormap jet
subplot(1,2,2); imagesc(atlas); axis off, grid on
saveas(gcf, [wfi_path, 'atlas_warped.fig'])

% Save the warped atlas
save([wfi_path, 'atlas_tform.mat'],'areanames','tform', 'invtform');
close all;