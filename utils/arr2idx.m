function idx = arr2idx(arr)

[M, N] = size(arr);
if M > 1 && N > 1, error('Array must be of size 1xN'); end
if M>N, arr = arr'; end


d = diff(arr);
start = find(d==1)+1;
stop = find(d==-1);


if start(1) > stop(1), start = [1 start]; end
if stop(end) < start(end), stop = [stop length(arr)]; end


assert(length(stop) == length(start), 'Must be same number of start and stop events')

idx = [start', stop'];
