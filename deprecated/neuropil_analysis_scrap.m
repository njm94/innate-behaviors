% plot neuropil correlations with grooming-responsive population

groom_pop = rastermap(pops(1,1):pops(1,2),:);

mean_groom = mean(groom_pop);

%% load 2p data

I = loadtiff('/home/user/teamshare/TM_Lab/nick/behavior/grooming/2p/ECL3_thy1/20240731/arm1/roi_1/bin2x2x1/data.tif');

%%
I2 = I;
for i = 1:length(cstat)
    I2(cstat{i}.ypix, cstat{i}.xpix,:) = nan;
end

%%

I3 = imresize(I2, 0.5, 'bilinear');

mean_groom_resamp = resample(groom_pop', size(I3,3), length(mean_groom))';


Ireshape = double(reshape(I3, size(I3,1)*size(I3,2), size(I3,3))');

corrmap = reshape(corr(mean_groom_resamp', Ireshape, 'rows', 'complete'), size(I3,1), size(I3,2));
figure, imagesc(corrmap)
colorbar
caxis([0 0.5])


%%

groom_pop_sub = mean_groom_resamp(:,6.4e3:6.75e3);
groom_pop_sub = groom_pop_sub - mean(groom_pop_sub,2);


[coeff,score,latent, ~,explained] = pca(groom_pop_sub);


%%


for i = 1%:length(long_eps)
    X = groom_pop(:,long_eps(i,1):long_eps(i,2));
    t = xt(X, fs, 2);
    for j = 1:size(X,1)
        avg_slope(j,:) = polyfit(t, X(j,:), 1);
    end
end
