% clear all
% close all
% clc

% generate minimum spanning tree
A = digraph(B);
% mst = minspantree(sparse(max(A(:)) - A),'method','kruskal');
% mst = minspantree(A);

% node properties
rngSz = [5,50];                     % range of node sizes
nodeSz = fcn_sz(sum(A,2),rngSz);    % node sizes

% visualization parameters
mapSz = 1000;    % size of map
sigma = 10;      % radius of smoothed colors

% generate thresholded network and union with mst
thr = 0.02;
A = double(threshold_proportional(A,thr) | mst | mst');

% make figure
[f,ax,nodeHandles,edgeHandles] = fcn_community_boundaries(A,lab,cmap,nodeSz,mapSz,sigma);