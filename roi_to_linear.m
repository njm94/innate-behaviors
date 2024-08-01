clear, clc

% addpath(genpath('C:\Users\user\Documents\Nick\grooming'));

[roi_file, roi_path] = uigetfile('*.tif','Select resonant scan ROI data.', 'Y:\nick\behavior\grooming\2p');
data_roi = loadtiff([roi_path, roi_file]);
frame_roi = uint16(mean(data_roi, 3));

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
[tform, moving_registered] = cpregister(imadjust(frame_roi), imadjust(frame_lin), tform_type);

figure, imshowpair(imadjust(frame_lin), imadjust(moving_registered))

% keep registration results in a separate variable since automated
% registration in the next step will overwrite the 'tform' variable
tform_estimate = tform; 


%% automated registration
% polish off registration with automated procedure using control point
% registration results as an initial transformation estimate
% visualize_registration(imadjust(frame_lin), [roi_path, 'ROIs.gif'], moving_registered)


warning('off', 'images:regmex:registrationOutBoundsTermination')
[optimizer,metric] = imregconfig("monomodal");
tform = imregtform(frame_roi, frame_lin, "affine", optimizer, metric, 'initialTrans', tform_estimate);
movingRegistered = imwarp(frame_roi,tform,'OutputView',imref2d(size(frame_lin)));
warning('on', 'images:regmex:registrationOutBoundsTermination')

figure, 
subplot(1,2,1), 
imshowpair(imadjust(frame_lin), imadjust(movingRegistered))
subplot(1,2,2), title('Final reg vs initial estimate')
imshowpair(imadjust(moving_registered), imadjust(movingRegistered))


%% Run this after accepting results of automated registration

moving_registered = movingRegistered;
figure, subplot(1,3,1), imshow(imadjust(frame_lin), []), 
subplot(1,3,2), imshow(moving_registered, [min(frame_roi(:)) prctile(frame_roi(:), 90)])
subplot(1,3,3), imshowpair(imadjust(frame_lin), moving_registered)
saveas(gcf, [roi_path, 'ROI_aligned.fig'])



% figure, imshowpair(imadjust(frame_lin), imadjust(moving_registered))
% savefig(gcf, [roi_path, 'ROI.fig'])
%% Save the roi to linear scan transformation matrix
save([roi_path, 'tform.mat'], 'tform', "frame_lin", "frame_roi", "moving_registered");
close all;



