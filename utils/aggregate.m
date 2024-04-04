function [bmat2, idx] = aggregate(bmat, W, fs)
% aggregate grooming events by merging all events that occur within W
% seconds of each other. Then return list of start and stop indexes

bmat2 = bmat>0;
if nargin < 3 || isempty(fs), fs = 90; end

if ~(nargin < 2 || isempty(W))
    W = round(W*fs); % convert seconds to samples
        
    bmat2 = movmax(bmat2, W);
    bmat2 = movmin(bmat2, W);
end




% figure, plot(bmat>0), hold on, plot(bmat2-1)

time_on = find(diff(bmat2)==1)+1;
time_off = find(diff(bmat2)==-1)+1;

% edge cases where grooming occurs at the beginning or end of trial
if bmat2(end), time_off = [time_off; length(bmat2)]; end
if bmat2(1), time_on = [1; time_on]; end

idx = [time_on; time_off]';
end