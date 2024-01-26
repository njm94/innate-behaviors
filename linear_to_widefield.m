clear, clc

addpath(genpath('Y:\nick\2p\code'));

[wfi_file, wfi_path] = uigetfile('*.tif','Select widefield reference image.', 'Y:\nick\behavior\grooming\2p')
frame_wfi = uint16(imflatfield(mean(data_wfi, 3), 100));

[lin_file, lin_path] = uigetfile('*.tif','Select 2-photon linear scan image.', 'Y:\nick\behavior\grooming\2p');
data_lin = loadtiff([lin_path, lin_file]);
frame_lin = uint16(mean(data_lin, 3));

%%

[tform, moving_registered] = cpregister(imadjust(frame_lin), imadjust(frame_wfi));

%%
gif_filename = [lin_path, 'ROIs.gif'];
visualize_registration(imadjust(frame_wfi), gif_filename, moving_registered)

%%
figure, subplot(1,2,1), imshow(imadjust(frame_wfi), []), 
subplot(1,2,2), imshowpair(imadjust(frame_wfi), moving_registered, 'blend')
saveas(gcf, [lin_path, 'linear_aligned.fig'])

%% Save the linear scan transformation matrix
save([lin_path, 'linearscan_tform.mat'], 'tform', "frame_wfi", "frame_lin", "moving_registered");
close all;

