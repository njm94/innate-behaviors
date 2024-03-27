function [f,ax,nodeHandles,edgeHandles] = fcn_community_boundaries(A,lab,cmap,nodeSz,mapSz,sigma)
% fcn_community_boundaries      makes a good looking network layout
%
%  [f,ax,nodeHandles,edgeHandles] = fcn_community_boundaries(A,lab,cmap,
%    nodeSz,mapSz,sigma) uses a force directed layout to plot locations of
%    nodes and edges. "A" is a connectivity matrix (must be sparse) with
%    community/cluster/category labels, "lab," and with cluster colors
%    specified by "cmap." Node sizes are encoded in the vector "nodeSz."
%    The variables "mapSz" and "sigma" are for plotting: "mapSz" is the
%    size of the background image (larger images take longer to generate --
%    maybe start with a number around 200 or 300 and go from there).
%    "sigma" is a smoothing parameter that determines how diffuse the
%    colors are behind each community.
%
%   Rick Betzel, Indiana University, 2019
%%
d = agreement(lab);             % agreement matrix for labels
coor = fcn_get_layout_coor(A);  % locations using matlab's graph layout function
xrng = ...                      % range of x coordinates (with buffer)
    [min(coor(:,1)) - 0.25*range(coor(:,1)),max(coor(:,1)) + 0.25*range(coor(:,1))];
yrng = ...                      % range of y coordinates (with buffer)
    [min(coor(:,2)) - 0.25*range(coor(:,2)),max(coor(:,2)) + 0.25*range(coor(:,2))];
%%
grid = zeros(mapSz,mapSz,max(lab)); % for storing node and edge locations
%%
for i = 1:max(lab)                  % loop over all categories
    gridtemp = zeros(mapSz);        % temporary grid
    idx = lab == i;                 % get current category
    ai = A(idx,idx);                % get subgraph
    coori = coor(idx,:);            % coordinates
    x = interp1(xrng,[1,mapSz],coori(:,1)); % map x coordinates to cells
    y = interp1(yrng,[1,mapSz],coori(:,2)); % same for y
    rx = round(x);                  % round to whole #
    ry = round(y);                  % ditto
    pts = (rx - 1)*mapSz + ry;
    gridtemp(pts) = 1;
    [u,v] = find(triu(ai));
    allpts = zeros(101*length(u),1);
    for iedge = 1:length(u)
        start = coori(u(iedge),:);
        finish = coori(v(iedge),:);
        xint = interp1([0,1],[start(1),finish(1)],linspace(0,1,101));
        yint = interp1([0,1],[start(2),finish(2)],linspace(0,1,101));
        x = interp1(xrng,[1,mapSz],xint);
        y = interp1(yrng,[1,mapSz],yint);
        rx = round(x);
        ry = round(y);
        edx = (1:101) + (iedge - 1)*101;
        allpts(edx) = (ry - 1)*mapSz + rx;
    end
    allpts = unique(allpts);
    gridtemp(allpts) = 1;
    grid(:,:,i) = gridtemp;
end
coor(:,1) = interp1(xrng,[1,mapSz],coor(:,1));
coor(:,2) = interp1(yrng,[1,mapSz],coor(:,2));
%%
sz = ceil(sqrt(mapSz));
x = linspace(-sz,sz);
y = linspace(-sz,sz);
[mx,my] = meshgrid(x,y);
kernel = exp(-(mx.^2 + my.^2)./(2*sigma.^2));
kernel = kernel/sum(kernel(:));
%%
tot = zeros(mapSz,mapSz,3);
colall = zeros(mapSz,mapSz,3,max(lab));
for i = 1:max(lab)
    gridsmooth = conv2(grid(:,:,i),kernel,'same');
    prob = gridsmooth/max(gridsmooth(:));
    tot(:,:,i) = prob;
    
    P = bsxfun(@minus,ones(mapSz,mapSz,3),prob);
    C = ones(mapSz,mapSz,3);
    for k = 1:3
        C(:,:,k) = cmap(i,k)*prob;
    end
    col = P + C;
    colall(:,:,:,i) = col;
end

%%
[~,idx] = max(tot,[],3);
cc = ones(mapSz,mapSz,3);
for i = 1:max(lab)
    mask = idx == i;
    [u,v] = find(mask);
    for j = 1:length(u)
        cc(u(j),v(j),:) = colall(u(j),v(j),:,i);
    end
end
%%
for i = 1:3
    cc(:,:,i) = conv2(cc(:,:,i),kernel,'same');
end
prct = round(mapSz*0.1);
cc(1:prct,:,:) = 1;
cc(:,1:prct,:) = 1;
cc(:,(end - prct + 1):end,:) = 1;
cc((end - prct + 1):end,:,:) = 1;
f = fcn_rickplot([2,2,5,5]);
ax = axes;
hold(ax,'on');
imagesc(flipud(rot90(cc)));
[ex,ey] = fcn_adjacency_plot(A & ~d,coor);
plot(ex,ey,'color',ones(1,3)*0.65);
nodeHandles = cell(max(lab),1);
edgeHandles = nodeHandles;
for i = 1:max(lab)
    idx = lab == i;
    ai = A(idx,idx);
    coori = coor(idx,:);
    [ex,ey] = fcn_adjacency_plot(ai,coori);
    edgeHandles{i} = plot(ex,ey,'color',cmap(i,:));
    nodeHandles{i} = scatter(coori(:,1),coori(:,2),nodeSz(idx),repmat(cmap(i,:),sum(idx),1),'o','filled');
    set(nodeHandles{i},'markeredgecolor','k');
end
axis image square;

function f = fcn_rickplot(pos)
f = figure(...
    'units','inches',...
    'position',pos,...
    'paperpositionmode','auto');

function [X,Y,Z] = fcn_adjacency_plot(aij,xy)
%FCN_ADJACENCY_PLOT     quick visualization tool
%
%   [X,Y,Z] = FCN_ADJACENCY_PLOT(AIJ,XY) takes adjacency matrix AIJ and
%   coordinates XY and generates output X,Y,Z that allows for quick
%   visualization of network using the PLOT command only. If no coordinates
%   are specified than the nodes are arranged on a ring.
%
%   Example:
%
%   >> load AIJ;                                % load your adjacency matrix
%   >> load COOR;                               % load [x,y,z] coordinates for each node
%   >> [x,y,z] = fcn_adjacency_plot(AIJ,COOR);  % call function
%   >> plot3(x,y,z);                            % plots network as a single line object
%
%   Richard Betzel, Indiana University, 2013

n = length(aij);
if nargin < 2
    xy = zeros(n,2);
    for i = 1:n
        xy(i,:) = [cos(2*pi*(i - 1)./n), sin(2*pi*(i - 1)./n)];
    end
end

if issymmetric(aij)
    aij = triu(aij,1);
end
[i,j] = find(aij);
[~, p] = sort(max(i,j));
i = i(p);
j = j(p);

X = [ xy(i,1) xy(j,1)]';
Y = [ xy(i,2) xy(j,2)]';
if size(xy,2) == 3
    Z = [ xy(i,3) xy(j,3)]';
end
if isfloat(xy) || nargout ~= 0
    X = [X; NaN(size(i))'];
    Y = [Y; NaN(size(i))'];
    if size(xy,2) == 3
        Z = [Z; NaN(size(i))'];
    end
end

X = X(:);
Y = Y(:);
if size(xy,2) == 3
    Z = Z(:);
end

function coor = fcn_get_layout_coor(adj)
[u,v,w] = find(triu(adj,1));
g = graph(u,v,w);
f = figure;
ph = plot(g,'layout','force','usegravity',true);
coor = [ph.XData',ph.YData'];
close(f);