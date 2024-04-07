clear, clc

data_root = 'Y:\nick\behavior\grooming\1p';
mice = {'ECR2_thy1', 'GER2_ai94', 'HYL3_tTA', 'IBL2_tTA'};



for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    dff_path = [data_root, filesep, mice{j}, filesep, 'outputs'];
    h = openfig([dff_path, filesep, getAllFiles(dff_path, '_dFF.fig')]);
    rightmove(:,:,j) = h.Children(8).Children.CData;
    leftmove(:,:,j) = h.Children(10).Children.CData;
    left(:,:,j) = h.Children(16).Children.CData;
    right(:,:,j) = h.Children(14).Children.CData;
    lick(:,:,j) = h.Children(12).Children.CData;
    elliptical(:,:,j) = h.Children(24).Children.CData;
    largeleft(:,:,j) = h.Children(22).Children.CData;
    largeright(:,:,j) = h.Children(20).Children.CData;
    dropright(:,:,j) = h.Children(2).Children.CData;
    bilateral(:,:,j) = h.Children(18).Children.CData;
    close(h)
end

%% overlay contours from diff mice
clc
nanmask = zeros(128, 128, length(mice));

for j = 1:length(mice)
    load([data_root, filesep, mice{j}, filesep, 'mask.mat'])
    nanmask(:,:,j) = mask;
end
nanmask(nanmask==0) = nan;





%%
clc, clear a v
tmp = [];

thresh = 80;
figure, axis off, hold on
for j = 1:length(mice)

    vars = ["lick", "right", "left", "elliptical", "largeright", "largeleft", "bilateral"];
         
    for i = 1:length(vars)    
        test = eval(vars(i));
        test = test .* nanmask;
        test = test(:,:,j);

        v = prctile(test(:), thresh);
        tmp = [tmp v];
        disp(v)
        a{i}(:,:,j) = test >= v;

        subplot(1,length(vars), i), axis off, hold on
            if v > 0
                contourf(flipud(test), [v v], 'FaceAlpha', 0.25)
    
                title(vars(i));
            end
%         end
    end
%     legend(labs, 'Location', 'Best')
end

%% compute pairwise dice
for i = 1:length(vars)
    count = 1;
    for j = 1:length(mice)
        for k = 1:length(mice)
            if k > j
                similarity(count, i) = dice(a{i}(:,:,j), a{i}(:,:,k));
                count = count + 1;
            end
        end
    end
end




[p,tbl,stats] = anova1(similarity);
[c,m,h,gnames] = multcompare(stats);

figure, boxplot(similarity, 'Colors', 'k'),
xticklabels(vars)
ylabel('Pairwise Dice Similarity Coefficient')
ax = gca;
ax.FontSize = 12;

title(['One-Way ANOVA, p = ', num2str(p, 2)])
