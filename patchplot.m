function patchplot(idx, yrange, color, alpha)
% Draw patches on plot timeseries plot
%
% idx       Nx2 array of start/end indices for n patches to draw
% yrange    ylimits (lower and upper edge of patches)
% color     color of patch
% alpha     opacity of patch

if nargin < 4, alpha = 0.5; end
x = zeros(4, size(idx,2));
for i = 1:size(idx,1)
    x(:,i) = [idx(i,1); idx(i,1); idx(i,2); idx(i,2)];
end
y = repmat([yrange(1); yrange(2); yrange(2); yrange(1);],[1, size(x,2)]);

patch(x, y, color, 'FaceAlpha', alpha, 'EdgeColor', 'None'), hold on

end

