clear, clc
%%
mice = {'Y:\nick\behavior\grooming\1p\ECR2_thy1', ...
    'Y:\nick\behavior\grooming\1p\GER2_ai94', ...
    'Y:\nick\behavior\grooming\1p\HYL3_tTA', ...
    'Y:\nick\behavior\grooming\1p\IBL2_tTA'};

load('allenDorsalMap.mat');

%%

[template_file, template_path] = uigetfile('*.tif','Select template image.', 'Y:\nick\behavior\grooming');
template = loadtiff([template_path, filesep, template_file]);

%%

tform = align_recording_to_allen(imadjust(template), {}, 1);

invT=pinv(tform.T); % invert the transformation matrix
invT(1,3)=0; invT(2,3)=0; invT(3,3)=1; % set 3rd dimension of rotation artificially to 0
invtform=tform; invtform.T=invT; % create the transformation with invtform as the transformation matrix
load('atlas.mat')


%%

imagereg = imwarp(imadjust(template),tform,'interp','nearest','OutputView',imref2d(size(dorsalMaps.dorsalMapScaled)));

% widefield image with atlas overlaid
figure
% imagesc(imflatfield(imagereg, 100)); colormap gray
imagesc(imagereg); colormap gray
axis equal off; hold on;
for p = 1:length(dorsalMaps.edgeOutline)
    plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1));
end
set(gca, 'YDir', 'reverse');
% saveas(gcf, [template_path, 'atlas_aligned.fig'])



%%

save([template_path, 'atlas_tform.mat'],'areanames','tform', 'invtform');
close all;

% %% Below uses the outputs from the regression to align the brain. 
% for i = 1:length(mice)
%     fit_is_not_good = true;
%     template = loadtiff([mice{i}, filesep, 'template.tif']);
%     load([mice{i}, filesep, 'mask.mat']);
%     % load regression figure to assist with registration of allen atlas
%     h = openfig([mice{i}, filesep, 'outputs', filesep, 'ridge_summary.fig']);
%     
%     % For now use position in h.Children to index variables. Not ideal, but
%     % couldn't figure out how to index by axis name
%     if contains(mice{i}, 'tTA')
%         vars = {'leftmove', 'rightmove', 'audio'};
%         audio = h.Children(24).Children.CData;
%     else
%         vars = {'leftmove', 'rightmove'};
%     end
%     rightmove = h.Children(2).Children.CData;
%     leftmove = h.Children(4).Children.CData;
% %     fullmat = h.Children(26).Children.CData;
% 
%     close(h);
% 
%     [fllx, flly] = get_loc_from_map(leftmove, 'centroid');
%     [flrx, flry] = get_loc_from_map(rightmove, 'centroid');
% %     figure
% %     for j = 1:length(vars)
% %         [x,y] = get_loc_from_map(eval(vars{j}), 'centroid');
% %         subplot(1,length(vars), j)
% %         imagesc(eval(vars{j})), hold on
% %         plot(x,y,'k*')
% %         title(vars{j})
% %     end
% %     
% 
%     % scale template and activity maps
% 
%     tmp = 65535 .* (template - min(template(:)))./ (max(template(:)) - min(template(:)));
%     leftmove = 65535 .* (leftmove - min(leftmove(:)))./ (max(leftmove(:)) - min(leftmove(:)));
%     rightmove = 65535 .* (rightmove - min(rightmove(:)))./ (max(rightmove(:)) - min(rightmove(:)));
%     img = cat(3, tmp.*uint16(mask), test.*mask, test2.*mask);
% 
% 
%     
%     while fit_is_not_good
%         tform = align_recording_to_allen(template, [], [], [fllx, flly; flrx, flry]);
%         figure('Position', [6 397 1268 538])
%     
%         tmp = {'leftmove', 'rightmove', 'template'};
%         
%         for ii = 1:3  
%             imagereg = imwarp(single(eval(tmp{ii})).*mask,tform,'interp','nearest','OutputView',imref2d(size(dorsalMaps.dorsalMapScaled)));
%             for p = 1:length(dorsalMaps.edgeOutline)
%                 for k = 1:size(dorsalMaps.edgeOutline{p})
%                     imagereg(round(dorsalMaps.edgeOutline{p}(k, 1)), round(dorsalMaps.edgeOutline{p}(k, 2))) = 0;
%                 end
%             end
%             subplot(1,3,ii)
%         
%             imagesc(single(imagereg))
%         end
%         answer = questdlg('Is the fit good?','','Yes','No','No');
%         switch answer
%             case 'Yes'
%                 fit_is_not_good = false;
%         end
%         close gcf
%     end
% 
%     sgdsg
% end
% 
% %%
% 
% 
% 
% %%
% 
% 
% test = ones(size(dorsalMaps.dorsalMapScaled));
% for p = 1:length(dorsalMaps.edgeOutline)
%     for k = 1:size(dorsalMaps.edgeOutline{p})
%         test(round(dorsalMaps.edgeOutline{p}(k, 1)), round(dorsalMaps.edgeOutline{p}(k, 2))) = nan;
%     end
% end
% % se = offsetstrel('disc',5,5);
% % test = imerode(test, se);
% test = imgaussfilt(test,0.1);
% 
% test = imwarp(test, invtform, 'interp', 'nearest', 'OutputView', imref2d(size(template)));
% 
% mask(mask==0) = nan;
% figure, pcolor(flipud(mask.*leftmove .* test)), shading flat
% %%
% clc
% h = openfig([mice{1}, filesep, 'outputs', filesep, 'ridge_summary.fig']);
% % h = openfig([mice{1}, filesep, 'outputs', filesep, '2024-03-08-07-56-29_dFF.fig']);
% 
% for i = 2:2:length(h.Children)
%     set(h.Children(i).Children, 'CData', h.Children(i).Children.CData.*test.*mask)
% end
% %%
% rightmove = h.Children(2).Children.CData;
% leftmove = h.Children(4).Children.CData;
% fullmat = h.Children(26).Children.CData;
% audio = h.Children(24).Children.CData;
% lick = h.Children(6).Children.CData;
% right = h.Children(8).Children.CData;
% left = h.Children(10).Children.CData;
% bilateral = h.Children(12).Children.CData;
% largeright = h.Children(14).Children.CData;
% largeleft = h.Children(16).Children.CData;
% elliptical = h.Children(18).Children.CData;
% 
% 
% figure
%         for ii = 1:3  
%             imagereg = imwarp(single(eval(tmp{ii})).*mask,tform,'interp','nearest','OutputView',imref2d(size(dorsalMaps.dorsalMapScaled)));
%             for p = 1:length(dorsalMaps.edgeOutline)
%                 for k = 1:size(dorsalMaps.edgeOutline{p})
%                     imagereg(round(dorsalMaps.edgeOutline{p}(k, 1)), round(dorsalMaps.edgeOutline{p}(k, 2))) = 1000;
%                 end
%             end
%             subplot(1,3,ii)
%         
%             imagesc(single(imagereg))
%         end
% 
% %%
% 
% 
% function [y,x] = get_loc_from_map(map, method, k)
% 
%     switch method
%         case 'max'
%             if nargin < 3 || isempty(k), k = 2; end
%             tmp = imgaussfilt(map, k);
%             [~,I] = max(tmp, [], 'all');
%             [x, y] = ind2sub(size(map), I);
%         case 'centroid'
%             if nargin < 3 || isempty(k), k = 99; end
%             P = prctile(map,k,"all");
%             [x,y] = ind2sub(size(map), find(map>P));
%             x = mean(x);
%             y = mean(y);
%             
%     end
% 
% end