function [bmat2, idx] = aggregate(bmat, W, fs)
% Aggregate binary events by merging all events that occur within W
% seconds of each other. Then return list of start and stop indices.
%
% Inputs:
%   bmat    (matrix, binary behavior matrix, [Time x Behavior])
%   W       (scalar, aggregation window size, in seconds)  
%   fs      (scalar, sampling rate of data)
%
% Outputs:
%   bmat2   (matrix, aggregated behavior matrix)
%   idx     (cell array, each cell is Nx2 array of start and stop index)
%
% Usage: 
%   [bmat2, idx] = aggregate(bmat, W, fs) 
%

% assume time along 1st dimension
if size(bmat,2) > size(bmat,1)
    bmat = bmat';
end



bmat2 = bmat>0;
if nargin < 3 || isempty(fs), fs = 90; end

if ~(nargin < 2 || isempty(W))
    W = round(W*fs); % convert seconds to samples
    % An even W ultimately results in a shifted index by 1 sample. Add 1 to
    % avoid this
    if ~mod(W,2)
        W = W+1;
    end
        
    bmat2 = movmax(bmat2, W);
    bmat2 = movmin(bmat2, W);
end

idx = arr2idx(bmat2);
end