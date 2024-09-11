function tform = sensory2atlas(image, pixelsize, sens, use_midline, tform_type)
%% Align recording to Allen CCF given an image 
% INPUTS
% image       - mean of the gcamp image
% scaling     - scaling of gcamp image
% sens        - sensory maps (HxWxN) where N is number of regions mapped
% use_midline - use line along midline for orientation (default true)
% tform_type  - default affine
%
% alignareas - additional areas in the form 'L area', with a space between
% L/R and area name. the user should be fairly certain of identifying these areas.
%
% affine - do an affine transformation or not. default: not an affine
% transformation. The alternative is a non-reflective similarity, aka
% rotation, translation, and scaling.
%
% OUTPUT
% tform - transformation details that can be run with imwarp for subsequent frames.

% Author: Shreya Saxena (2018)
%%
if nargin<2 || isempty(pixelsize), pixelsize=1000*8.2/1024; end
% if nargin<3 || isempty(sens), ; end
if nargin<4 || isempty(use_midline), use_midline = true; end
if nargin<5 || isempty(tform_type), tform_type = 'similarity'; end

load('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat');
load('C:\Users\user\Documents\Nick\grooming\utils\atlas.mat')
         
list = {'R SSp-ul','L SSp-ul', ...
    'R SSp-ll', 'L SSp-ll', ...
    'R SSp-bfd', 'L SSp-bfd', ...
    'R AUD', 'L AUD'};


% apply scaling tform
scalingFactor = pixelsize./dorsalMaps.desiredPixelSize;

% Create the isotropic scaling matrix
S = [scalingFactor 0 0;
     0 scalingFactor 0;
     0 0 1];

Stform = affine2d(S);

% get image in coordinate size from the start
sImage = imwarp(image, Stform, 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
for p = 1:size(sens,3)
    sSens(:,:,p) = imwarp(sens(:,:,p), Stform, 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
end


if use_midline
    handle = figure;
    imagesc(sImage); axis equal off;hold on; colormap gray
    title('Click on points along midline from top to bottom, then press enter','fontsize',12);
    [x,y] = getline(handle);

    ang = -atan((x(2)-x(1))./abs(y(2)-y(1)));
    
    % Create the 2D rotation matrix
    R = [cos(ang) -sin(ang) 0;
         sin(ang)  cos(ang) 0;
         0 0 1];
    Rtform = affine2d(R);
    
    rsImage = imwarp(sImage, Rtform, 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
    rsSens = imwarp(sSens, Rtform, 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));
%     rsImage = imrotate(sImage, rad2deg(-ang));
%     rsSens = imrotate(sSens, rad2deg(-ang));

    % Rotate the midline coordinates to get the vertical line bregma is on
    rMidline = R(1:2, 1:2) * [x, y]';
    bregma_x = mean(rMidline(1,:));

    % Bregma center for dorsalMapScaled from Shreya Saxena code
    bregref_x= 293.9059; 
    close(handle)
end


figure
for i = 1:size(sens,3)
    imshowpair(double(rsImage), double(rsSens(:,:,i)));
    [x,y] = ginput(1);

    [indx,~] = listdlg('PromptString','Select the region',...
    'SelectionMode','single','ListString',list);

    alignareas{i} = list{indx};
    points(i, 1) = x;
    points(i, 2) = y;
    close(gcf)
end


points_refs=zeros(length(alignareas),2);
for i=1:length(alignareas)
    sidestr=alignareas{i}(1);
    try
        areaid=dorsalMaps.labelTable{strcmp(dorsalMaps.labelTable.abbreviation,alignareas{i}(3:end)),'id'}+1;

    catch
        error('Please name additional areas in the form L/R area, with a space between L/R and area name');
    end
    dimsmap=size(dorsalMaps.dorsalMap);
    if areaid == 99
        % If areaid is equal to 99, then this is auditory cortex. Combine
        % auditory subregions for this one
        auditory_areas = [100, 114, 107];
        x = [];
        y = [];
        for jj = 1:length(auditory_areas)
            [xt,yt]= find(dorsalMaps.dorsalMapScaled==auditory_areas(jj));
            x = [x; xt];
            y = [y; yt];
        end
    else
        [x,y]=find(dorsalMaps.dorsalMapScaled==areaid);
    end
    
    if strcmp(sidestr,'R')
        x=x(y>=dimsmap(2)/2); y=y(y>=dimsmap(2)/2);
    elseif strcmp(sidestr,'L')
        x=x(y<dimsmap(2)/2); y=y(y<dimsmap(2)/2);
    end

    if areaid == 99
        % If auditory cortex is the reference, use the most medial points
        if strcmp(sidestr,'R')
            points_refs(i,:)=[mean(y), min(x)];
        elseif strcmp(sidestr,'L')
            points_refs(i,:) = [mean(y), max(x)];
        end
    else
        points_refs(i,:)=[mean(y) mean(x)];
    end
end

if use_midline
%     xShift = mean([mean(points_refs(:,1) - points(:,1)), bregref_x - bregma_x]);
%     xShift = mean([points_refs(:,1) - points(:,1); bregref_x - bregma_x]);
%     xShift = bregref_x - bregma_x;
    xShift = mean(points_refs(:,1) - points(:,1));
    yShift = mean(points_refs(:,2) - points(:,2));


%     figure, imagesc(rsImage), colormap gray
%     axis equal off; hold on;
%     hold on
% 
%     for p = 1:length(dorsalMaps.edgeOutline)
%         plot(dorsalMaps.edgeOutline{p}(:, 2)-xShift, dorsalMaps.edgeOutline{p}(:, 1)-yShift, 'w', 'LineWidth', 2);
%     end
% 
%     for p = 1:size(points,1)
%         plot(points(p, 1), points(p,2), 'rx')
%     end

    % create translation matrix
    T = [1 0 xShift;
        0 1 yShift;
        0 0 1];

    % combine all the transformation matrices into a composite affine
    % transformation matrix

    % change the direction of rotation since imwarp rotates in the
    % opposite direction for some reason    
    R(2,1) = -R(2,1);
    R(1,2) = -R(1,2);
    TRS = T*R*S;
    tform = affinetform2d(TRS);
end

imagereg = imwarp(image, tform, 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));

figure
imagesc(imagereg); colormap gray, hold on

for p = 1:length(dorsalMaps.edgeOutline)
    plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'w', 'LineWidth', 2);
end


end