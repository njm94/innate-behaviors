b = 6;

figure, imagesc(behavior_frames{b}{1}(:,:,1))
[x,y] = ginput(1);
y = round(y);
x=round(x);
close(gcf)

sizes = cellfun(@(x) size(x,3), behavior_frames{b});
[~, idx] = sort(sizes);

test = nan(length(sizes), max(sizes));
for i = 1:length(sizes)
    test(i,1:sizes(idx(i))) = squeeze(behavior_frames{b}{idx(i)}(y,x,:));
end

figure, imagesc(test), 
colormap(redblue)
colorbar
caxis([-4 4])



