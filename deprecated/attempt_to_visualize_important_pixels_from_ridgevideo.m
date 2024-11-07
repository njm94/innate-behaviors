N = size(Ume,1);
for i = 1:round(sqrt(N))
    if ~rem(N,i)
        disp(i)
    end

end

%%

[coeff,score,latent] = pca(UU);

%%

figure, plot(latent)

%%
