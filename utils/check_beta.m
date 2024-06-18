function cBeta = check_beta(cVar, regLabels, regIdx, U, dimBeta, Vc, mask, visual)
%   required inputs
% cVar:             variable of interest. Must be included in 'regLabels'
% regLabels:        regressor labels
% regIdx:           regressor index
% U:                spatial components from SVD
% dimBeta:          outputs from regression
% Vc:               temporal components from SVD
%   optional inputs
% mask:             brain mask
% visual:           visualize using the GUI

if nargin < 7 || isempty(mask), mask = ones(128, 128); end
if nargin < 8 || isempty(visual), visual = 'True'; end

% find beta weights for current variable
cIdx = regIdx == find(ismember(regLabels,cVar));
U = reshape(U, [], size(Vc,1)); 
cBeta = U * dimBeta(cIdx, :)';
cBeta = reshape(cBeta, size(mask,1), size(mask,2), []); 

cBeta = cBeta .* mask;
if visual, compareMovie(cBeta); end

end