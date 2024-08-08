clear, clc

load('C:\Users\user\Documents\Nick\grooming\utils\atlas.mat')
load('C:\Users\user\Documents\Nick\grooming\utils\allen_map\allenDorsalMap.mat');

startpath = 'Y:\nick\behavior\grooming\2p';
[r_tform, roipath] = uigetfile('*.mat','Select roi transformation.', startpath);
roi_tform = load([roipath, r_tform]);

[suite2p_file, suite2p_path] = uigetfile('*.mat', 'Select Suite2p output file.', roipath);
load([suite2p_path, suite2p_file])

[l_tform, startpath] = uigetfile('*.mat','Select linear scan transformation.', roipath);
lscan_tform = load([startpath, l_tform]);

[w_tform, startpath] = uigetfile('*.mat','Select widefield transformation.', startpath);
wfi_tform = load([startpath, w_tform]);

%% Unpack frames from loaded files

frame_roi = roi_tform.frame_roi;
frame_lin = lscan_tform.frame_lin;
frame_wfi = lscan_tform.frame_wfi;
imreg = imwarp(frame_wfi, wfi_tform.tform,'interp','nearest','OutputView',imref2d(size(dorsalMaps.dorsalMapScaled)));
% imreg = imwarp(frame_wfi, wfi_tform.tform);

 
%%

show_figs = true;
if show_figs
    h1 = figure(1); imshow(imadjust(frame_roi)), axis off; hold on
    h2 = figure(2); imshow(imadjust(frame_lin)), axis equal off; hold on
    h3 = figure(3); imshow(imadjust(frame_wfi)), axis equal off; hold on
    h4 = figure(4); imshow(imadjust(imreg), []), axis equal off; hold on
%     h4 = figure(4); imshow(imadjust(frame_wfi), []), axis equal off; hold on

end



area_values = struct2array(areanames);
area_names = fieldnames(areanames);

xShift = wfi_tform.tform.A(1,3);
yShift = wfi_tform.tform.A(2,3);

for i = 1:size(F,1)
    if iscell(i,1)
        x0 = double(stat{i}.med(2));
        y0 = double(stat{i}.med(1));        
        [x1, y1] = transformPointsForward(roi_tform.tform, x0, y0);        
        [x2, y2] = transformPointsForward(lscan_tform.tform, x1, y1);        
        [x3, y3] = transformPointsForward(wfi_tform.tform, x2, y2);

        % some issue with transformation not being applied properly to
        % dorsalMaps variable. As a workaround, translate the dorsalMaps
        % with the final translation and check those area names
%         dMap = dorsalMaps.dorsalMapScaled;
% 
%         xShift = wfi_tform.tform.A(1,3);
%         yShift = wfi_tform.tform.A(2,3);
%         dMap = imtranslate(dMap, [-xShift, -yShift], 'nearest');
        n_loc(i) = area_names(area_values==dorsalMaps.dorsalMapScaled(round(y3), round(x3)));
    

        if show_figs
            figure(1), plot(x0, y0, 'ro')
            figure(2), plot(x1, y1, 'r.', 'MarkerSize', 10)
            figure(3), plot(x2, y2, 'r.', 'MarkerSize', 10)
            figure(4), plot(x3, y3, 'r.', 'MarkerSize', 10)
        end
    else
        n_loc(i) = {'NULL'};
    end
end
if show_figs
    figure(4)
    for p = 1:length(dorsalMaps.edgeOutline)
        plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'w', 'LineWidth', 1.75);
    end
    set(gca, 'YDir', 'reverse');
    saveas(gcf, [roipath, 'atlas_aligned.fig'])
end

%%
save([suite2p_path, suite2p_file], 'n_loc', '-append')
disp("Done")


