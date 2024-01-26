clear, clc

% addpath(genpath('C:\Users\user\Documents\Nick\grooming'));

[roi_file, roi_path] = uigetfile('*.tif','Select resonant scan ROI data.', 'Y:\nick\behavior\grooming\2p');
data_roi = loadtiff([roi_path, roi_file]);
frame_roi = uint16(mean(data_roi, 3));

[lin_file, lin_path] = uigetfile('*.tif','Select 2-photon linear scan image.', 'Y:\nick\behavior\grooming\2p');
data_lin = loadtiff([lin_path, lin_file]);
frame_lin = uint16(mean(data_lin, 3));

%%

[tform, moving_registered] = cpregister(imadjust(frame_roi), imadjust(frame_lin));
%%

visualize_registration(imadjust(frame_lin), [roi_path, 'ROIs.gif'], moving_registered)
figure, subplot(1,2,1), imshow(frame_lin, []), 
subplot(1,2,2), imshowpair(frame_lin, moving_registered, 'blend')

%% Save the roi to linear scan transformation matrix
save([roi_path, 'tform.mat'], 'tform', "frame_lin", "frame_roi", "moving_registered");
close all;



