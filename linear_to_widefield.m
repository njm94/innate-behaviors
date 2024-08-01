clear, clc

addpath(genpath('Y:\nick\2p\code'));

[wfi_file, wfi_path] = uigetfile('*.tif','Select widefield reference image.', 'Y:\nick\behavior\grooming\2p');
data_wfi = loadtiff([wfi_path, wfi_file]);
frame_wfi = uint16(imflatfield(mean(data_wfi, 3), 100));

[lin_file, lin_path] = uigetfile('*.tif','Select 2-photon linear scan image.', 'Y:\nick\behavior\grooming\2p');
data_lin = loadtiff([lin_path, lin_file]);
frame_lin = uint16(mean(data_lin, 3));

%% do control point registration
% tform_type can be one of the following
% 'nonreflectivesimilarity'
% 'similarity'
% 'affine'
% 'projective'
% 'pwl'

tform_type = 'affine';
[tform, moving_registered] = cpregister(imadjust(frame_lin), imadjust(frame_wfi), tform_type);
figure, imshowpair(imadjust(frame_wfi), imadjust(moving_registered))

% keep registration results in a separate variable since automated
% registration in the next step will overwrite the 'tform' variable
tform_estimate = tform; 

%% automated registration
% polish off registration with automated procedure using control point
% registration results as an initial transformation estimate
% 

warning('off', 'images:regmex:registrationOutBoundsTermination')
[optimizer,metric] = imregconfig("multimodal");
tform = imregtform(frame_lin, frame_wfi, "affine", optimizer, metric, 'initialTrans', tform_estimate);
movingRegistered = imwarp(frame_lin,tform,'OutputView',imref2d(size(frame_wfi)));
warning('on', 'images:regmex:registrationOutBoundsTermination')

figure, 
subplot(1,2,1), title('Final reg vs WFI')
imshowpair(imadjust(frame_wfi), imadjust(movingRegistered))
subplot(1,2,2), title('Final reg vs initial estimate')
imshowpair(imadjust(moving_registered), imadjust(movingRegistered))

%% Run this after accepting results of automated registration

moving_registered = movingRegistered;
figure, subplot(1,3,1), imshow(imadjust(frame_wfi), []), 
subplot(1,3,2), imshow(imadjust(moving_registered), [])
subplot(1,3,3), imshowpair(imadjust(frame_wfi), imadjust(moving_registered))
saveas(gcf, [lin_path, 'linear_aligned.fig'])

%% Save the linear scan transformation matrix
save([lin_path, 'linearscan_tform.mat'], 'tform', "frame_wfi", "frame_lin", "moving_registered");
close all;

