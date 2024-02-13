function arr = idx2arr(idx, L)

arr = zeros(1, L);
for i = 1:size(idx, 1)
    arr(idx(i,1):min([idx(i,2), L])) = 1;
end
