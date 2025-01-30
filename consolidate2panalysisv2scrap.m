clc, clear, close all
fileID = fopen('expt1_datalist.txt','r');
addpath('C:\Users\user\Documents\Nick\grooming\utils')

formatSpec = '%s';
% mp_list1 = {'Y:\nick\behavior\grooming\2p\ETR2_thy1\20231113143925'; ...
%     'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231113155903'; ...
%     'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231115174148';
%     };

mp_list = {%'Y:\nick\behavior\grooming\2p\ETR2_thy1\20231113143925';
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231113155903';
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231115174148';
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240731';
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240802';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240729';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240731';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240802';
    'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240729';
    'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240802'};

% mp_list = {'Y:\nick\behavior\grooming\2p\ETR2_thy1\20231113143925'; ...
%     'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231113155903'; ...
%     'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231115174148'; ...
%     'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240729'; ...
%     'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240731';
%     'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240802';
%     'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240729';
%     'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240731';
%     'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240802';
%     'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240729';
%     'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240802'};

mp_list = fix_path(mp_list);
data_list = mp_list;
current_mouse = '';
cluster_data = fix_path('Y:\nick\behavior\grooming\20241114092737_behavior_clustering.mat');
dmetric = 'cosine';

addpath(genpath('C:\Users\user\Documents\Nick\grooming\utils'))
try load('C:\Users\user\Documents\Nick\grooming\utils\allen_map\allenDorsalMap.mat');
catch
    load('/home/user/Documents/grooming/utils/allen_map/allenDorsalMap.mat');
end

fs = 90;
include_boris = true;
%%

% Plot the Allen Atlas
figure(1), axis equal off; hold on; set(gca, 'YDir', 'reverse');

for p = 1:length(dorsalMaps.edgeOutline)-2 % -2 to ignore olfactory bulbs
    plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k', 'LineWidth', 2);
end
%%

all_possible_labs = {'FLR', 'FLL', 'Left', 'Right', 'Elliptical', 'Left Asymmetric', 'Right Asymmetric', 'Elliptical Right', 'Elliptical Left'};
consolidated_labs = {'FL', 'Unilateral', 'Elliptical', 'Asymmetric', 'Ellip-Asymm'};
total_neurons = 0;
for i = 1:length(mp_list)
    
    mouse_root = fileparts(mp_list{i});
    
    if ~strcmp(mouse_root, current_mouse)
        load([mouse_root, filesep, 'dalsa', filesep, 'atlas_tform.mat'])
    end

    boris_file = [mp_list{i}, filesep, getAllFiles(mp_list{i}, 'events.tsv')];



    if ~isfile([mp_list{i}, filesep, 'Nresample.mat'])
        load([mp_list{i}, filesep, getAllFiles(mp_list{i}, 'Fclean.mat')]);
        
        [events, b_idx, ~, vid_end, cluster_labels] = get_labels_from_clustering_results(cluster_data, boris_file, include_boris);

        try Nresample = resample(N, vid_end, size(N, 2), 'Dimension', 2);
        catch
            disp('Data too big. Splitting into halves and resampling each half separately')
            Nresample = resamplee(N', size(events,1), size(N,2))';
        end
    
        save([mp_list{i}, filesep, 'Nresample.mat'], 'Nresample', 'nloc', 'fs', 'cstat', 'tforms', '-v7.3')
    else
        disp('Loading resampled neuron data')
        load([mp_list{i}, filesep, 'Nresample.mat'])
    end


%%

    load([mp_list{i}, filesep, 'Nresample.mat'])
    clear x3 y3
    for j = 1:length(cstat)
        x0 = double(cstat{j}.med(2));
        y0 = double(cstat{j}.med(1));        
        [x1, y1] = transformPointsForward(tforms{cstat{j}.use_tform}.roi_to_linear.tform, x0, y0);        
        [x2, y2] = transformPointsForward(tforms{cstat{j}.use_tform}.linear_to_wfield.tform, x1, y1);        
        [x3(j), y3(j)] = transformPointsForward(tforms{cstat{j}.use_tform}.wfield_to_atlas.tform, x2, y2);

    end
    total_neurons = total_neurons + length(x3);
    figure(1), scatter(x3, y3, 'k.')



    %%
    % Load BORIS file
%     [events, b_idx, ~] = read_boris([mp_list{i}, filesep, getAllFiles(mp_list{i}, 'events.tsv')]);
    [events, b_idx, ~, vid_end, cluster_labels] = get_labels_from_clustering_results(cluster_data, boris_file, include_boris);

    % load DLC tracks
    vel = readmatrix([mp_list{i}, filesep, getAllFiles(mp_list{i}, '_vel.csv')]);
    vel = vel(1:vid_end,:);
    flrv = sum(vel(:,4:5).^2, 2).^0.5;
    fllv = sum(vel(:,7:8).^2, 2).^0.5;
    flrthresh = flrv>mean(flrv) + std(flrv);
    fllthresh = fllv>mean(fllv) + std(fllv);
    
    % Create the behavior matrix for hierarchical clustering
    clear Bmean event_table labs tmp_Bmean
    count = 1;
    labs = {};
    for j = 1:size(events,2)
        % Get rid of events non-motor events
        if contains(events.Properties.VariableNames{j}, 'Drop') || contains(events.Properties.VariableNames{j}, 'Video') || contains(events.Properties.VariableNames{j}, 'Lick')
            continue
        end
        event_table(:, count) = table2array(events(:,j));
        Bmean(:,count) = mean(Nresample(:,logical(event_table(:,count))),2);
        labs{count} = events.Properties.VariableNames{j};
        count = count + 1;        
    end
    
    % Remove grooming movements from forelimb movements, then add to the
    % behavior matrix
    flrthresh(any(event_table, 2)) = 0;
    fllthresh(any(event_table, 2)) = 0;
    
    Bmean = cat(2, Bmean, mean(Nresample(:, logical(flrthresh)),2));
    Bmean = cat(2, Bmean, mean(Nresample(:, logical(fllthresh)),2));
    labs = [labs, 'FLR', 'FLL'];

    %%
    figure, hold on
    for j = 1:size(event_table,2)
        tmp = corr(event_table(:,j), Nresample');
        plot(sort(tmp))
    end
    legend(labs, 'Location', 'Best')


    t = xt(Nresample, fs, 2);
    figure, hold on
    for j = 1:size(event_table,2)
%         subplot(size(event_table,2),1,j)
        [tmp, I] = sort(corr(event_table(:,j), Nresample'));
        plot(t,zscore(Nresample(I(end),:))-j*5)
    end

    %% 
    eth_table = array2table(event_table, 'VariableNames', labs(1:size(event_table,2)));
    plot_ethogram(eth_table, labs(1:size(event_table,2)), fs)


    %%

    % Do the hierarchical clustering
    
    Z = linkage(Bmean', 'average', dmetric);

    figure, subplot(3,2,1)
    [H, T, outperm] = dendrogram(Z,'Labels', labs);  % Plot the dendrogram
    ylabel(['Distance (',dmetric,')'])
    xticks([])
        
    % sort neurons by anatomical location
    label_column = ones(size(Bmean,1), 2);
    anat = sort(unique(nloc));
    I = [];
    uni_idx = strcmp(labs, 'FLR') | strcmp(labs, 'FLL');

    for j = 1:length(anat)
        Itmp = find(strcmp(nloc,anat{j}));
        [~, Itmp_sorted] = sort(mean(Bmean(Itmp,uni_idx),2));
        I = [I; Itmp(Itmp_sorted)];
    end
    
    anat_parent = unique(cellfun(@(x) x(1:3), anat, 'UniformOutput', false));
    for j = 1:length(anat_parent)
        label_column(contains(nloc(I), anat_parent{j}), 1) = j;  
    end
    
    for j = 1:length(anat)
        label_column(contains(nloc(I), anat{j}), 2) = j;    
    end
    
    % label_column
    
    subplot(3,2,[3,5])
    imagesc(Bmean(I, outperm)),
    xticks(1:length(labs))
    xticklabels(labs(outperm))
    % c=colorbar;
    % c.Label.String = 'Z-score';
    caxis([-1 3])
    colormap(bluewhitered())
    ylabel('Neuron')
    freezeColors
    % 
    subplot(3,2,[4, 6]), imagesc(label_column), colormap default, axis off


    % Not all possible behaviors are observed in each session. To compare
    % across sessions, map the relationships between neuronal population
    % distances across behaviors to a matrix which contains all possible
    % behaviors. If some behaviors are not present, fill those with NaNs

    for j=1:length(all_possible_labs)
        if any(strcmp(all_possible_labs{j}, labs))
            tmp_Bmean(j,:) = Bmean(:, strcmp(all_possible_labs{j}, labs));
        else
            tmp_Bmean(j,:) = nan(1, size(Bmean,1));
        end
    end

    sim_matrix(:,:,i) = 1-squareform(pdist(tmp_Bmean, dmetric));



end


%%


dat = mean(sim_matrix,3,'omitnan');


[trueM,Q]=community_louvain(dat);
count = 1;
clear M ari
for i = 1:size(sim_matrix,3)
    if isnan(mean(sim_matrix(:,:,i), 'all'))
        continue
    end
    M(:, count) = community_louvain(sim_matrix(:,:,i), [], [], 'negative_asym');
    ari(count) = rand_index(trueM, M(:, count));
    count = count + 1;
end

figure, imagesc(dat), 
caxis([0.5 1]),
colormap(flipud(colormap(gray)))
pt = find(diff(trueM));

line([0.5 0.5], [0 pt+0.5], 'Color', [0 0.4470 0.7410], 'LineWidth', 5)
line([pt+0.5 pt+0.5], [0 pt+0.5], 'Color', [0 0.4470 0.7410], 'LineWidth', 5)
line([0.5 pt+0.5], [0.5 0.5], 'Color', [0 0.4470 0.7410], 'LineWidth', 5)
line([0.5 pt+0.5], [pt+0.5 pt+0.5], 'Color', [0 0.4470 0.7410], 'LineWidth', 5)

line([pt+0.5 pt+0.5], [pt+0.5 9.5], 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 5)
line([9.5 9.5], [pt+0.5 9.5], 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 5)
line([pt+0.5 9.5], [9.5 9.5], 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 5)
line([pt+0.5 9.5], [pt+0.5 pt+0.5], 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 5)

axis equal off
% exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'behavior_cluster_2p.png']), 'Resolution', 300)
%%
figure, axis off; caxis([0.5 1]), colormap(flipud(colormap(gray))); colorbar
% exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'behavior_cluster_2p_colorbar.png']), 'Resolution', 300)