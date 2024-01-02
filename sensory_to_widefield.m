clear, clc

addpath(genpath('Y:\nick\2p\code'));

[wfi_file, wfi_path] = uigetfile('*.tif','Select widefield reference image.', 'Y:\nick\2p');
data_wfi = loadtiff([wfi_path, wfi_file]);
frame_wfi = uint16(imflatfield(mean(data_wfi, 3), 100));
%%

[sens_file, sens_path] = uigetfile('*.tif','Select sensory stimulus image.', 'Y:\nick\2p');
data_sens = loadtiff([sens_path, sens_file]);
frame_sens = uint16(mean(data_sens, 3));

[dff_file, dff_path] = uigetfile('*.raw','Select averaged sensory dFF.', 'Y:\nick\2p');
dff_sens = open_raw([dff_path, dff_file], [], [], [128 128]);

%% transform frames (acquisition is upside down)

frame_wfi = imrotate(frame_wfi, 180);
%%
frame_sens = imrotate(frame_sens, 180);
dff_sens = imrotate(dff_sens, 180);


%% image registration
fixed = frame_wfi;
moving = frame_sens;

[tform, moving_registered] = cpregister(imadjust(frame_sens), imadjust(frame_wfi));
figure, imshowpair(frame_wfi, moving_registered)
savefig([sens_path,'tform_results.fig'])
save([sens_path,'tform.mat'], 'tform')

%% get the response map

stim_time = round(size(dff_sens, 3)/2);


for p = 1:4;
    sensory_response = mean(dff_sens(:,:,stim_time:stim_time+9), 3);
    sensory_response(sensory_response<0) = 0;
    sensory_response = sensory_response .^ p; % visually separate response from background
    sensory_map = imwarp(sensory_response, tform, 'OutputView',imref2d(size(frame_wfi)));
    
    figure, imshowpair(frame_wfi, sensory_map)
    
    sensory_map = uint16(max(double(frame_wfi(:))) .* (sensory_map-min(sensory_map(:))) ./ (max(sensory_map(:)) - min(sensory_map(:))));
    imwrite(sensory_map, [sens_path, 'response_', num2str(p),'.tiff'])
end

%%



