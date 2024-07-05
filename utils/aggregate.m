function [bmat2, idx] = aggregate(bmat, W, fs)
% aggregate grooming events by merging all events that occur within W
% seconds of each other. Then return list of start and stop indexes

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