% sensory_to_atlas
clear, clc, close all
[comp_file, comp_path] = uigetfile('*.tif','Select sensory composite image.', 'Y:\nick\behavior\grooming\2p');
data_comp = loadtiff([comp_path, comp_file]);
load('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat');


%%


frame_wfi = data_comp(:,:,1);
sensory_maps = data_comp(:,:,2:end);

FOV = 8.2; % mm
scaling = 1000 .* FOV/size(frame_wfi,1);
use_midline = true;


tform = sensory2atlas(frame_wfi, scaling, sensory_maps, use_midline);


%%


list = {'R SSp-ul','L SSp-ul', ...
    'R SSp-ll', 'L SSp-ll', ...
    'R SSp-bfd', 'L SSp-bfd'};




figure
for i = 2:size(data_comp,3)
    imshowpair(data_comp(:,:,1), data_comp(:,:,i));
    [x,y] = ginput(1);
%     x = round(x); y = round(y);

    [indx,~] = listdlg('PromptString','Select the region',...
    'SelectionMode','single','ListString',list);

    regions{i-1} = list{indx};
    region_coords(i-1, 1) = x;
    region_coords(i-1, 2) = y;
end

%%
% frame_wfi = uint16(data_comp(:,:,1));
% sensory_maps = double(data_comp(:,:,2:end));
% 
% 
% 
% % % get max location for reference
% % for i = 1:size(sensory_maps, 3)
% %     tmp = imgaussfilt(sensory_maps(:,:,i), 10);
% %     bw_sensory(:,:,i) = regiongrowing(tmp);
% %     bw_centroids(i,:) = round(regionprops(bw_sensory(:,:,i)).Centroid);
% % end
% 
% 
% %%
% sensory_maps = double(sum(sensory_maps, 3)) .* double(max(frame_wfi(:)));
% 
% figure, imshowpair(frame_wfi, sensory_maps)
% 
% % clear test;
% test = cat(3,double(frame_wfi), sensory_maps);
% 
% %%
% figure, imagesc(frame_wfi), colormap gray, hold on
% imagesc(double(bw_sensory(:,:,1)) .* double(max(frame_wfi(:))))
% %%



%% Do registration and transformation on atlas

% tform = align_recording_to_allen_v2(frame_wfi, regions, 1, region_coords); % align <-- input any function of data here
invT=pinv(tform.T); % invert the transformation matrix
invT(1,3)=0; invT(2,3)=0; invT(3,3)=1; % set 3rd dimension of rotation artificially to 0
invtform=tform; invtform.T=invT; % create the transformation with invtform as the transformation matrix
load('atlas.mat')
% atlas=imwarp(atlas,invtform, 'interp','nearest', 'OutputView', imref2d(size(frame_wfi))); % setting the 'OutputView' is important
% atlas=round(atlas);


%%



imagereg = imwarp(frame_wfi,tform,'interp','nearest','OutputView',imref2d(size(dorsalMaps.dorsalMapScaled)));
sensoryreg = imwarp(sensory_maps, tform, 'interp', 'nearest', 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
% widefield image with atlas overlaid
figure
% imagesc(imflatfield(imagereg, 100)); colormap gray
% imagesc(imagereg); colormap gray
for i = 1:size(sensory_maps,3)
    subplot(1, size(sensoryreg,3), i)
    imshowpair(imagereg, sensoryreg(:,:,i))
    axis equal off; hold on;
    for p = 1:length(dorsalMaps.edgeOutline)
        plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'w', 'LineWidt', 2);
    end
    set(gca, 'YDir', 'reverse');
end
%%
saveas(gcf, [comp_path, 'atlas_aligned.fig'])

%%

% 
save([comp_path, 'atlas_tform.mat'],'areanames','tform', 'invtform');
% close all;


