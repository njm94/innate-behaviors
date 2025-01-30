clc, clear, close all
fileID = fopen('expt1_datalist.txt','r');
addpath('C:\Users\user\Documents\Nick\grooming\utils')

formatSpec = '%s';
% mp_list1 = {'Y:\nick\behavior\grooming\2p\ETR2_thy1\20231113143925'; ...
%     'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231113155903'; ...
%     'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231115174148';
%     };



mp_list = {'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240731';
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240802';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240729';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240731';
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231113155903';
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231115174148';};
mp_list = fix_path(mp_list);
current_mouse = '';
cluster_data = fix_path('Y:\nick\behavior\grooming\20241114092737_behavior_clustering.mat');

dmetric = 'cosine';
fs = 90;
aggregation_sz = 3; % window of time to aggregate behaviors
addpath(genpath('/home/user/Documents/grooming/utils'))
try load('C:\Users\user\Documents\Nick\grooming\utils\allen_map\allenDorsalMap.mat');
catch
    load('/home/user/Documents/grooming/utils/allen_map/allenDorsalMap.mat');
end

include_boris = true;
%%

all_possible_labs = {'FLR', 'FLL', 'Left', 'Right', 'Elliptical', 'Left Asymmetric', 'Right Asymmetric', 'Elliptical Right', 'Elliptical Left'};
consolidated_labs = {'FL', 'Unilateral', 'Elliptical', 'Asymmetric', 'Ellip-Asymm'};

region_locs = {};
region_neurons = [];
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
    uiopen([mouse_root, filesep, 'dalsa', filesep, 'atlas_aligned.fig'],1), hold on
    scatter(x3, y3,50, '.', 'MarkerEdgeColor', [0 200 150]/255)

    %%
    % Load BORIS file
    [events, b_idx, ~, vid_end, cluster_labels] = get_labels_from_clustering_results(cluster_data, boris_file, include_boris);

    % load DLC tracks
    vel = readmatrix([mp_list{i}, filesep, getAllFiles(mp_list{i}, '_vel.csv')]);
    vel = vel(1:vid_end,:);
    flrv = sum(vel(:,4:5).^2, 2).^0.5;
    fllv = sum(vel(:,7:8).^2, 2).^0.5;
    flrthresh = flrv>mean(flrv) + 3.5*std(flrv);
    fllthresh = fllv>mean(fllv) + 3.5*std(fllv);
    
    % Create the behavior matrix for hierarchical clustering
    clear Bmean event_table labs tmp_Bmean
    count = 1;
    labs = {};
    for j = 1:size(events,2)
        % Get rid of events non-motor events
        if contains(events.Properties.VariableNames{j}, 'Drop') || contains(events.Properties.VariableNames{j}, 'Video') %|| contains(events.Properties.VariableNames{j}, 'Lick')
            continue
        end
        event_table(:, count) = table2array(events(:,j));
        Bmean(:,count) = mean(Nresample(:,logical(event_table(:,count))),2);
        labs{count} = events.Properties.VariableNames{j};
        count = count + 1;        
    end
    
    % aggregate drop events into a single column
    idx = contains(events.Properties.VariableNames, 'Drop');
    Drop = any(table2array(events(:,idx)),2);
    events = removevars(events, idx);
    events = addvars(events, Drop);


    % get index of long grooming events
    bmat = any(table2array(removevars(events, 'Drop')), 2);
    [episodes, idx] = aggregate(bmat, aggregation_sz);
    episode_durations = diff(idx, 1, 2)/90;
    long_eps = [];
    for jj = 1:length(episode_durations)
        if episode_durations(jj) >= 10
            long_eps = cat(1, long_eps, idx(jj,:));
        end
    end

    % Remove grooming movements from forelimb movements, then add to the
    % behavior matrix
    flrthresh(aggregate(any(event_table, 2), 3, fs)) = 0;
    fllthresh(aggregate(any(event_table, 2),3, fs)) = 0;


    
    Bmean = cat(2, Bmean, mean(Nresample(:, logical(flrthresh)),2));
    Bmean = cat(2, Bmean, mean(Nresample(:, logical(fllthresh)),2));
    labs = [labs, 'FLR', 'FLL'];

    %%

    load([mp_list{i}, filesep, 'rastermap_order.mat'])
    
    h2 = figure('Position', [180 251 1458 536]); hold on
    rastermap = Nresample(isort+1,:);
    
    imagesc(imadjust(rastermap)), colormap(flipud(colormap(gray)))
    % offset = size(Nresample,1);
    Move = flrthresh|fllthresh;
    states = {'Start', 'Right', 'Left', 'Elliptical', ...
        'Right Asymmetric', 'Left Asymmetric', ...
        'Elliptical Right', 'Elliptical Left', 'Stop', 'Lick', 'Drop', 'Move'};
    plot_ethogram(addvars(events, Move), states, 1, size(Nresample,1), 50)
    axis tight
    
    axis off
 
    if i == 1
        pops = [540 663;
            1 40;
            52 115;
            178 216;
            ];
        zoom_x = [45725 88684];
    elseif i == 2
        pops = [1 32;
            70 99];
        
    elseif i == 3
        pops = [300 380;
            385 419;
            ];
            % 1 30];
    elseif i == 4
        pops = [144 199;
            28 71];
    elseif i == 5
        pops = [198 237;
            1 105];
    elseif i == 6
        pops = [115 137;
            1 66];
    end


            cols = [0 0 1;
            1 0 0;
            1 0 1;
            0 1 1];

    nloc_reordered = nloc(isort+1);
    for jj = 1:size(rastermap,1)
        
        if any(contains(region_locs, nloc_reordered(jj)))
            region_idx = find(contains(region_locs, nloc_reordered(jj)), 1);
        else
            region_locs = [region_locs, nloc_reordered(jj)];
            disp(length(region_locs))
            region_idx = length(region_locs);
            region_neurons = cat(1, region_neurons, zeros(1, 3));
        end
        

        if any(pops(1,1):pops(1,2) == jj)
            region_neurons(region_idx, 1) = region_neurons(region_idx,1) + 1;
        elseif any(pops(2,1):pops(2,2) == jj)
            region_neurons(region_idx, 2) = region_neurons(region_idx,2) + 1;
        else
            region_neurons(region_idx,3) = region_neurons(region_idx,3) + 1;
        end
    end

    fs = 90;
    rastermap(:,[1:fs, end-fs:end]) = nan; % get rid of the ringing artifact in the start and end due to upsampling
    % move_neurons = mean(rastermap(pops(2,1):pops(2,2),:));
    % for jj = 1:size(rastermap,1)
    %     move_corr(jj) = corr(rastermap(i,:), move_neurons);
    % end

    for jj = 1:size(pops,1)
        line([0 0], pops(jj,:), 'Color', cols(jj,:), 'LineWidth', 5)
        plot(30*mean(rastermap(pops(jj,1):pops(jj,2),:))- 100*jj, 'LineWidth', 1, 'Color', cols(jj,:))
    end
    line([1 fs*60], [-500 -500], 'Color', [0 0 0], 'LineWidth', 2)
    % 
    % vline(zoom_x, 'k-')
    % exportgraphics(h2, fix_path(['Y:\nick\behavior\grooming\figures\thy1_mouse_ethogram_raster.png']), 'Resolution', 300)
    
    % axis([zoom_x, ylim])
    % line([zoom_x(1)+1 zoom_x(1)+1+fs*20], [-500 -500], 'Color', [0 0 0], 'LineWidth', 2)
    % axis([64253    75989 ylim])
    % line([64253 64253 + 90*10], [-50 -50], 'Color', [0 0 0], 'LineWidth', 2)
    % exportgraphics(h2, fix_path(['Y:\nick\behavior\grooming\figures\thy1_ethogram_raster_zoom.png']), 'Resolution', 300)
    
    uiopen([mouse_root, filesep, 'dalsa', filesep, 'atlas_aligned.fig'],1), hold on
    scatter(x3, y3,25, '.', 'MarkerEdgeColor', [0 0 0])
    for jj = 1:size(pops,1)
        idx = isort+1;
        scatter(x3(idx(pops(jj,1):pops(jj,2))), y3(idx(pops(jj,1):pops(jj,2))), 16, 'MarkerFaceColor', cols(jj,:), 'MarkerEdgeColor', 'k')
    end
    mean_groom_pop = mean(rastermap(pops(1,1):pops(1,2),:));
    mean_move_pop = mean(rastermap(pops(2,1):pops(2,2),:));
    for jj = 1:size(long_eps,1)      
        try
            long_groom{i}{jj} = mean_groom_pop(long_eps(jj,1)-5*fs:long_eps(jj,2)); 
            long_move{i}{jj} = mean_move_pop(long_eps(jj,1)-5*fs:long_eps(jj,2)); 
        catch
            continue
        end
    end


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

long_groom_mat = [];
ep_length = [];
for i = 1:length(long_groom)
    for j = 1:length(long_groom{i})
        tmp_length = size(long_groom{i}{j},2);
        ep_length = [ep_length tmp_length];
        if isempty(long_groom_mat)
            long_groom_mat = long_groom{i}{j};
        elseif tmp_length > size(long_groom_mat,2)
            size_dif = tmp_length - size(long_groom_mat,2);
            long_groom_mat = [long_groom_mat, nan(size(long_groom_mat,1), size_dif)];
            long_groom_mat = [long_groom_mat; long_groom{i}{j}];
        elseif tmp_length < size(long_groom_mat,2)
            size_dif = size(long_groom_mat,2) - tmp_length;
            long_groom_mat = [long_groom_mat; 
                [long_groom{i}{j} nan(1, size_dif)]];
        else
            long_groom_mat = [long_groom_mat; long_groom{i}{j}];
        end
    end
end



long_move_mat = [];
ep_length = [];
for i = 1:length(long_move)
    for j = 1:length(long_move{i})
        tmp_length = size(long_move{i}{j},2);
        ep_length = [ep_length tmp_length];
        if isempty(long_move_mat)
            long_move_mat = long_move{i}{j};
        elseif tmp_length > size(long_move_mat,2)
            size_dif = tmp_length - size(long_move_mat,2);
            long_move_mat = [long_move_mat, nan(size(long_move_mat,1), size_dif)];
            long_move_mat = [long_move_mat; long_move{i}{j}];
        elseif tmp_length < size(long_move_mat,2)
            size_dif = size(long_move_mat,2) - tmp_length;
            long_move_mat = [long_move_mat; 
                [long_move{i}{j} nan(1, size_dif)]];
        else
            long_move_mat = [long_move_mat; long_move{i}{j}];
        end
    end
end





%%
figure
[~,I] = sort(ep_length, 'descend');
fs = 90;
t = xt([long_groom_mat(1,:), 0], fs, 2)-5;
num_ep = size(long_groom_mat,1):-1:1;

% pcolor crops off the last row and column. add a column to the end to fix
% this. There is one row of nans already (due to one grooming episode that
% occurs at the beginning of the trial) so no need to add a row
pcolor(t, num_ep, [long_groom_mat(I,:), zeros(size(I))']), 

shading flat
caxis([-1 3])
colormap(bluewhitered())


p50 = round(prctile(ep_length(ep_length>0), 50));

line(([p50 p50]/fs)-5, ylim, 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '--')
line([0 0], ylim, 'Color', [0 0 0], 'LineWidth', 2)
line([51 60], [3 3], 'Color', [0 0 0], 'LineWidth', 2)
xticklabels([])
yticklabels([])
xticks([])
yticks([])
set(gca, 'YDir', 'reverse')

% exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'longgroom_heatmap_2p.png']), 'Resolution', 300)


%%
figure, axis off, colorbar, caxis([-1 3]), colormap(bluewhitered())
% exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'longgroom_heatmap_2p_colorbar.png']), 'Resolution', 300)
%%

figure, hold on
shadedErrorBar(t(1:p50), mean(long_groom_mat(:,1:p50), 'omitnan'), std(long_groom_mat(:,1:p50), 'omitnan'), 'k', 1)
line([0 0], ylim, 'Color', [0 0 0], 'LineWidth', 2)
axis([-5 30 -1 2.5])

% xticklabels([])
% yticklabels([])
% exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'longgroom_average2p.png']), 'Resolution', 300)

%%

% remove the empty row now
empty_idx = sum(long_groom_mat,2,'omitnan') == 0;
long_groom_mat(empty_idx,:) = [];
ep_length(empty_idx) = [];
t = xt(long_groom_mat, fs, 2)-5;
t_pre = t<0;
% t_on = abs(t)<=1;
t_early = t>=0 & t<5;
t_late = t>=5;

dff_pre = mean(long_groom_mat(:, t_pre), 2, 'omitnan');
% dff_on = mean(long_groom_mat(:, t_on), 2, 'omitnan');
dff_early = mean(long_groom_mat(:,t_early), 2, 'omitnan');
for i = 1:length(ep_length)
    dff_late(i) = mean(long_groom_mat(i, ep_length(i)-(5*fs):end), 2, 'omitnan');
end


% plot_data = [dff_pre, dff_on, dff_early, dff_late];
plot_data = [dff_pre, dff_early, dff_late'];
figure, boxplot(plot_data, 'Colors', 'k', 'Symbol', '')
hold on
swarmchart(repmat([1 2 3], size(dff_early,1), 1), plot_data, 'k', 'XJitterWidth', 0.25)
ylabel('Mean \DeltaF/F_0 (\sigma)')
xticklabels({'Pre 5s', 'First 5s', 'Last 5s'})
ax = gca;
ax.FontSize = 14;


% [p,~,stats] = anova1([dff_on, dff_early, dff_late])
% [c,m,h,gnames] = multcompare(stats);

clc
% Create a table with the data
tbl = array2table(plot_data, 'VariableNames', {'Pre', 'Early', 'Late'});

% Define the repeated measures model
rm = fitrm(tbl, 'Pre-Late ~ 1', 'WithinDesign', [1 2 3]'); % 3 time periods

% Run repeated-measures ANOVA
ranova_results = ranova(rm);

% Display the results
disp(ranova_results);

% Perform pairwise comparisons
pairwise_results = multcompare(rm, 'Time', 'ComparisonType', 'bonferroni');
disp(pairwise_results);

% saveas(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'binned_states2p.svg']))

%%
blen = 5;

clear p
for i = 1:size(long_groom_mat,1)
    
    X = long_groom_mat(i,round(blen*fs):ep_length(i));
    t = xt(X, fs);
    avg_slope(i,:) = polyfit(t, X, 1);
    
    % figure, plot(t, X)
    % f = fit(t,X,'exp1')
    % hold on, plot(t, f)
    % hfdhd
    % linslope(i) = p(1);
    % disp(p)
end
% close all
figure, boxplot(avg_slope(:,1), 'Colors', 'k', 'Symbol', '')
hold on, hline(0, 'k--')
swarmchart(ones(size(avg_slope,1),1), avg_slope(:,1), 'k', 'XJitterWidth', 0.5)
ylabel('Slope of best-fit line')
xticks([])
ax = gca; 
ax.FontSize = 14;


[h, p, ci, stats] = ttest(avg_slope(:,1), 0, 'Tail', 'left');

% Display results
disp(['p-value: ', num2str(p)]);
disp(['Test statistic (t): ', num2str(stats.tstat)]);
disp(['95% confidence interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);

if h == 1
    disp('The mean is significantly less than zero.');
else
    disp('The mean is not significantly less than zero.');
end



% saveas(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'average_slope2p.svg']))

%%

% figure
parent_region = cellfun(@(x) x(1:end-2), region_locs, 'UniformOutput', false);
unique_region = unique(parent_region);
clear consolidated_region_neurons
for i = 1:length(unique_region)
    idx = contains(region_locs, unique_region{i});
    consolidated_region_neurons(i,:) = sum(region_neurons(idx,:), 1);
end

groom_prop = zeros(size(dorsalMaps.dorsalMapScaled));
move_prop = zeros(size(dorsalMaps.dorsalMapScaled));
none_prop = zeros(size(dorsalMaps.dorsalMapScaled));
mp = round(size(dorsalMaps.dorsalMapScaled,2)./2);
for j = 1:length(unique_region)
    % If the value is negative, it is on the right side
    addmask = dorsalMaps.dorsalMapScaled == abs(areanames.([unique_region{j},'_L']));
    addmask(:,mp:end) = 0;
%     if areanames.(region_locs{j}) < 0
%         addmask(:, 1:mp) = 0;
%     else
%         addmask(:,mp:end) = 0;
%     end
    groom_prop = groom_prop + addmask.* (consolidated_region_neurons(j,1) ./ sum(consolidated_region_neurons(j,:)));
    move_prop = move_prop + addmask.* (consolidated_region_neurons(j,2) ./ sum(consolidated_region_neurons(j,:)));
    none_prop = none_prop + addmask.* (consolidated_region_neurons(j,3) ./ sum(consolidated_region_neurons(j,:)));
end

cols = [0 0 1; 
    1 0 0;
    0 0 0];


% figure, 
% h = bar((consolidated_region_neurons ./ sum(consolidated_region_neurons,2)), 'stacked');
% % Apply colors
% title(sum(consolidated_region_neurons,2))
% for i = 1:length(h)
%     h(i).FaceColor = 'flat';       % Enable flat coloring
%     h(i).CData = repmat(cols(i, :), size(consolidated_region_neurons, 1), 1); % Assign colors
% end
% legend({'Groom + Move', 'Move only', 'Neither'})
% ylabel('Proportion')
% xticklabels(unique_region)
% % xlabel('Episode duration (s)')
% % ax = gcf;

figure,  hold on
for i = 1%:3
    % subplot(1,3,i), hold on
    switch i
        case 1, imagesc(-groom_prop), caxis([-0.5 0]); colormap(bluewhitered()),  freezeColors
        case 2, imagesc(move_prop), colormap(bluewhitered()), caxis([0 0.5]); freezeColors
        case 3, imagesc(none_prop), colormap(flipud(colormap(gray))), caxis([0 1]), %h3=colorbar;
    end
    for p = 1:58%length(dorsalMaps.edgeOutline)
        % if p == 60, continue; end
        plot(dorsalMaps.edgeOutline{p}(:,2), dorsalMaps.edgeOutline{p}(:,1), 'k', 'LineWidth', 2)
    end
            axis equal off
        set(gca, 'YDir', 'reverse');
end

%%

exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'grooming_proportion2p.png']), 'Resolution', 300)

%%
