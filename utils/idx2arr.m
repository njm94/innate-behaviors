function arr = idx2arr(idx, L)
if nargin<2 || isempty(L), L=max(idx(:)); end

arr = zeros(1, L);
for i = 1:size(idx, 1)
    arr(idx(i,1):min([idx(i,2), L])) = 1;
end
