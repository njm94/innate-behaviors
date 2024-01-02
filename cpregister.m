function [tform, movingRegistered] = cpregister(moving, fixed, tformtype)
% Generic rigid registration pipeline
%
% Input:       
%           fixed, moving       (images for registration)
%           tformtype           (Default 'NonreflectiveSimilarity')
% Output:       
%           tform               (transformation matrix)
%           movingRegistered    (registered Image)
%

if nargin<3 || isempty(tformtype)
    tformtype = 'NonreflectiveSimilarity';
end

movingTmp = double(moving);
fixed = double(fixed);
movingTmp = (movingTmp-min(movingTmp(:))) ./ (max(movingTmp(:)) - min(movingTmp(:)));
fixed = (fixed-min(fixed(:))) ./ (max(fixed(:)) - min(fixed(:)));

[mp,fp] = cpselect(movingTmp,fixed,'Wait',true);
tform = fitgeotrans(mp,fp,tformtype);
movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));

end