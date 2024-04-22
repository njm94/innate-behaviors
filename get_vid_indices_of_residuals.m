


v = VideoReader('Y:\nick\behavior\grooming\1p\ECR2_thy1\20231114151923\ECR2_thy1_20231114151923_trim.mp4');

%%
loidx = find(lores);
loframes = [];
for i = 1:length(loidx)
    loframes = cat(3, loframes, read(v, loidx(i)));
end

%%
hiidx = find(hires);
hiframes = [];
for i = 1:length(hiidx)
    hiframes = cat(3, hiframes, read(v, hiidx(i)));
end

%%