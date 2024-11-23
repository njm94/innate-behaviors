% dice on individual mice


clc, clear

addpath('/home/user/Documents/grooming/utils')
addpath('C:\Users\user\Documents\Nick\grooming\utils')

fileID = fopen('expt1_datalist.txt','r');

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
data_list{1} = fix_path(data_list{1});
current_mouse = '';

fs = 90 ;

thy1_idx = 1:7;
ai94_idx = 8:13;
hyl3_idx = 14:19;
ibl2_idx = 20:25;

spon_index = false(length(data_list{1}),1);
spon_index([1,2,8,9,14,15,20,21]) = true;

save_average_across_days = true;

load(fix_path('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat'));



%%
mice = {'thy1_idx', 'ai94_idx', 'hyl3_idx', 'ibl2_idx'};
figure(1)
figure(2)
for ii = 1:length(mice)
    clear right* left* ellip* large* lick

    for j = eval(mice{ii})
        % ignore spontaneous trials
        if spon_index(j), continue; end
    
    
        mouse_root = fileparts(data_list{1}{j});
        load([mouse_root, filesep, 'mask.mat'])
        load([mouse_root, filesep, 'atlas_tform.mat'])
    
        dff_path = [data_list{1}{j}, filesep, 'outputs'];
        dff_fig = getAllFiles(dff_path, '_dFF.fig');
    
        % If there are multiple versions, use the most recent one
        if size(dff_fig,1) > 1
            dff_fig = sort(dff_fig);
            dff_fig = dff_fig{end};
        end
        h = openfig([dff_path, filesep, dff_fig]);
        for i = 1:length(h.Children)
            switch(h.Children(i).Title.String)
                case 'RightMove' 
                    rightmove(:,:,j) = h.Children(i).Children(end).CData;
                case 'LeftMove' 
                    leftmove(:,:,j) = h.Children(i).Children(end).CData;
                case 'largebilateral'
                    bilateral(:,:,j) = h.Children(i).Children(end).CData;
                case 'Elliptical'
                    elliptical(:,:,j) = h.Children(i).Children(end).CData;
                case 'Elliptical Right'
                    ellip_right(:,:,j) = h.Children(i).Children(end).CData;
                case 'Elliptical Left'
                    ellip_left(:,:,j) = h.Children(i).Children(end).CData;
                case 'Left Asymmetric'
                    largeleft(:,:,j) = h.Children(i).Children(end).CData;
                case 'Right Asymmetric'
                    largeright(:,:,j) = h.Children(i).Children(end).CData;
                case 'Left'
                    left(:,:,j) = h.Children(i).Children(end).CData;
                case 'Right'
                    right(:,:,j) = h.Children(i).Children(end).CData;
                case 'Lick'
                    lick(:,:,j) = h.Children(i).Children(end).CData;
                otherwise
                    continue
            end
        end
    
        close(h)
    end
    
    %%
    
    right(right==0) = nan;
    left(left==0) = nan;
    ellip_left(ellip_left==0) = nan;
    ellip_right(ellip_right==0) = nan;
    largeright(largeright==0) = nan;
    largeleft(largeleft ==0) = nan;
    elliptical(elliptical==0)=nan;
    lick(lick==0)=nan;
    
    
    %%
    clear avg_signal binary_map labs
    % close all
    vars = ["lick", "right", "left", "elliptical", "largeright", "largeleft", ...
        "ellip_right", "ellip_left"];
    thresh = 90;
    labcount = 1;
    
    figure(1)
    for i = 1:length(vars)
        behavior_var = eval(vars(i));
        bcount = 1;
        for j = 1:size(behavior_var,3)
            tmp = behavior_var(:,:,j);
            % if no behaviors in the session, continue
            if isnan(mean(tmp, [1 2], 'omitnan')), continue; end
    
            avg_signal{i}(bcount, 1) = mean(tmp, [1 2], 'omitnan');
    
            v = prctile(tmp(:), thresh);
            binary_map{i}(:,:,bcount) = tmp >= v;
            data_map{i}(:,:,bcount) = tmp;
            labs{labcount} = vars(i);
            labcount = labcount + 1;
            bcount = bcount+1;
        end
        subplot(length(mice), length(vars), length(vars)*(ii-1)+i),
        imagesc(mean(binary_map{i},3))
        clim([0 1]), colormap(bluewhitered())
    
        axis off, hold on, set(gca, 'YDir', 'reverse');
        title(vars(i));
        for p = 1:length(dorsalMaps.edgeOutline)
            plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k', 'LineWidth', 1);
            xticks([])
            yticks([])
        end
    end
    
    
    
    %% compute pairwise dice
    test = catcell(3, binary_map);
    dicemat = zeros(size(test,3), size(test,3));
    
    for i = 1:size(test,3)
        for j = 1:size(test,3)
            dicemat(j,i) = dice(test(:,:,j), test(:,:,i));
        end
    end
    %%
    
    figure(2), 
    subplot(length(mice), 1, ii)
    imagesc(dicemat), hold on
    colormap(flipud(colormap(gray)))
    bcount = 0.5;
    for i = 1:length(vars)
        bcount = bcount + size(binary_map{i},3);
    
        hline(bcount, 'r-')
        vline(bcount, 'r-')
    end
    xticks(cumsum(cellfun(@(x) size(x, 3), binary_map))-2)
    yticks(cumsum(cellfun(@(x) size(x, 3), binary_map))-2)

    xticklabels(vars)
    yticklabels(vars)
    colormap(flipud(colormap(gray)))
    c = colorbar;
    c.Label.String = 'Dice Similarity Coefficient';
    title('DFF binary')
end

%%

figure, 
imagesc(dicemat), hold on
for i = 1:length(vars)-1
    hline(5*(i)+0.5, 'r-')
    vline(5*(i)+0.5, 'r-')
end
xticks(2.5:5:size(dicemat,1))
yticks(2.5:5:size(dicemat,1))
xticklabels(vars)
yticklabels(vars)
% caxis([0 1]), 
% colormap gray
colormap(flipud(colormap(gray)))
c = colorbar;
c.Label.String = 'Dice Similarity Coefficient';
title('DFF binary')