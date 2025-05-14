function timeseries = getTimeseries(data, coords, radius)
% Get ROI timeseries from image stack. If no ROIs are specified return 
% global median
%
% Inputs:
%   data        (required, image stack with frames along the 3rd dimension)
%   coords      (required, ROI coordinates, N x 2 [x y])
%   radius      (optional, px surrounding coordinate location, default 2)


if nargin == 1
    timeseries = squeeze(nanmedian(nanmedian(data,2),1));
else
    if nargin<3
        radius = 2;
    end
    
    T = size(data, 3);
    R = size(coords, 1);
    timeseries = zeros(T, R);
    
    for i = 1:R
        regiony  = round(coords(i,2)-radius):round(coords(i,2)+radius);
        regionx  = round(coords(i,1)-radius):round(coords(i,1)+radius);
        timeseries(:, i) = squeeze(nanmean(nanmean(...
                    data(regiony,regionx,:), ...
                    2),1));

    end
    


end

end