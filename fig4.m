% Figure 4

clc, clear, close all
fileID = fopen('expt1_datalist.txt','r');
addpath('C:\Users\user\Documents\Nick\grooming\utils')
addpath(genpath('/home/user/Documents/grooming/utils'))
formatSpec = '%s';



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


mp_list = fix_path(mp_list);
current_mouse = '';
cluster_data = fix_path('Y:\nick\behavior\grooming\20241114092737_behavior_clustering.mat');

fs = 90;
aggregation_sz = 3; % window of time to aggregate behaviors

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
ncorr = [];
gcorr = [];
ncorrm = [];
gcorrm = [];

plot_atlas = false;
plot_rasters = false;
skip_region_data = false;

%%

for i = 1:length(mp_list)
    %% Loading
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



    %% Plot the Allen Atlas
if plot_atlas
    figure(1), axis equal off; hold on; set(gca, 'YDir', 'reverse');
    
    for p = 1:length(dorsalMaps.edgeOutline)-2 % -2 to ignore olfactory bulbs
        plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k', 'LineWidth', 2);
    end
end

    %% Load neuron data and rastermap embeddings
    load([mp_list{i}, filesep, 'Nresample.mat'])
    clear x3 y3
    for j = 1:length(cstat)
        x0 = double(cstat{j}.med(2));
        y0 = double(cstat{j}.med(1));        
        [x1, y1] = transformPointsForward(tforms{cstat{j}.use_tform}.roi_to_linear.tform, x0, y0);        
        [x2, y2] = transformPointsForward(tforms{cstat{j}.use_tform}.linear_to_wfield.tform, x1, y1);        
        [x3(j), y3(j)] = transformPointsForward(tforms{cstat{j}.use_tform}.wfield_to_atlas.tform, x2, y2);
    end
    fs = 90;
    total_neurons = total_neurons + length(x3);

    load([mp_list{i}, filesep, 'rastermap_order.mat'])    
    rastermap = Nresample(isort+1,:);
    nloc_reordered = nloc(isort+1);

    % get rid of the ringing artifact in the start and end due to upsampling
    rastermap(:,[1:fs, end-fs:end]) = nan; 
    % z-score again after getting rid of those artifacts
    rastermap = (rastermap - nanmean(rastermap,2))./nanstd(rastermap, [], 2);

    %% Index rastermap for neuronal populations
    % these are the neuronal populations manually identifed from the 
    % rastermap embeddings. The first row of pops are the indices of the
    % grooming population. Other highlighted populations are for
    % visualization in the example figure but most aren't used for anything
    if contains(mp_list{i}, 'ECL3') 
        if contains(mp_list{i}, '0731')
            gpop = 540:663;

            % this is the example so more data for plots
            pops = [540 663;
                1 40;
                52 115;
                178 216;
                ];
            zoom_x = [45725 88684];
            plot_rasters = true;
        elseif contains(mp_list{i}, '0802')
            gpop = 1:32;
        else
            skip_region_data = true;
        end
    elseif contains(mp_list{i}, 'IDR3') 
        if contains(mp_list{i}, '0729')
            gpop = 300:380;
        elseif contains(mp_list{i}, '0731')
            gpop = 144:199;
        else
            skip_region_data = true;
        end
    elseif contains(mp_list{i}, 'ETR3')
        if contains(mp_list{i}, '5903')
            gpop = 198:237;
        elseif contains(mp_list{i}, '4148')
            gpop = 115:137;
        else
            skip_region_data = true;
        end
    elseif contains(mp_list{i}, 'RR3')
        skip_region_data = true;
        if contains(mp_list{i}, '0802')
            % This mouse has massive ringing artifact on one of the microscope
            % arms. Exclude those neurons
            rastermap(1:140,:) = [];
            gpop = [359:384, 799:834];
        else
            continue
        end
    end
    %% Load behavior data
    [events, b_idx, ~, vid_end, cluster_labels] = get_labels_from_clustering_results(cluster_data, boris_file, include_boris);

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

    %% Load DLC data and get neuronal means for each behavior
    
    % load DLC tracks
    vel = readmatrix([mp_list{i}, filesep, getAllFiles(mp_list{i}, '_vel.csv')]);
    vel = vel(1:vid_end,:);
    flrv = sum(vel(:,4:5).^2, 2).^0.5;
    fllv = sum(vel(:,7:8).^2, 2).^0.5;
    flrthresh = flrv>mean(flrv) + 1*std(flrv);
    fllthresh = fllv>mean(fllv) + 1*std(fllv);
    
    % Remove grooming movements from forelimb movements, then add to the
    % behavior matrix
    flrthresh(aggregate(any(event_table, 2), 3, fs)) = 0;
    fllthresh(aggregate(any(event_table, 2),3, fs)) = 0;

    Bmean = cat(2, Bmean, mean(Nresample(:, logical(flrthresh)),2));
    Bmean = cat(2, Bmean, mean(Nresample(:, logical(fllthresh)),2));
    labs = [labs, 'FLR', 'FLL'];

    %% Plot rasters for figure
    
    if plot_rasters
        cols = [0 0 1;
            1 0 0;
            1 0 1;
            0 1 1];
        h2 = figure('Position', [180 251 1458 536]); hold on
        imagesc(imadjust(rastermap)), colormap(flipud(colormap(gray)))
        % offset = size(Nresample,1);
        Move = flrthresh|fllthresh;
        states = {'Start', 'Right', 'Left', 'Elliptical', ...
            'Right Asymmetric', 'Left Asymmetric', ...
            'Elliptical Right', 'Elliptical Left', 'Stop', 'Lick', 'Drop', 'Move'};
        plot_ethogram(addvars(events, Move), states, 1, size(Nresample,1), 50)
        axis tight
        for jj = 1:size(pops,1)
        line([0 0], pops(jj,:), 'Color', cols(jj,:), 'LineWidth', 5)
        plot(30*mean(rastermap(pops(jj,1):pops(jj,2),:))- 100*jj, 'LineWidth', 1, 'Color', cols(jj,:))
        end
        line([1 fs*60], [-500 -500], 'Color', [0 0 0], 'LineWidth', 2)
        plot_rasters = false;

        % vline(zoom_x, 'k-')
        % exportgraphics(h2, fix_path(['Y:\nick\behavior\grooming\figures\thy1_mouse_ethogram_raster.png']), 'Resolution', 300)
        
        % axis([zoom_x, ylim])
        % line([zoom_x(1)+1 zoom_x(1)+1+fs*20], [-500 -500], 'Color', [0 0 0], 'LineWidth', 2)
        % axis([64253    75989 ylim])
        % line([64253 64253 + 90*10], [-50 -50], 'Color', [0 0 0], 'LineWidth', 2)
        % exportgraphics(h2, fix_path(['Y:\nick\behavior\grooming\figures\thy1_ethogram_raster_zoom.png']), 'Resolution', 300)
        % axis off
    end
    
    %% Get locations of grooming populations
    
    if skip_region_data % No region information or no long groom episodes
        continue
    else
        skip_region_data = false;

        % get locations of grooming populations
        for jj = 1:size(rastermap,1)
            
            if any(contains(region_locs, nloc_reordered(jj)))
                region_idx = find(contains(region_locs, nloc_reordered(jj)), 1);
            else
                region_locs = [region_locs, nloc_reordered(jj)];
                region_idx = length(region_locs);
                region_neurons = cat(1, region_neurons, zeros(1, 2));
            end
            
            % if there is a grooming neuron in region X, add to first
            % column. otherwise, add to second
            
            if any(gpop == jj)
                region_neurons(region_idx, 1) = region_neurons(region_idx, 1) + 1;
            else
                region_neurons(region_idx, 2) = region_neurons(region_idx, 2) + 1;
            end
        end
       
        % compute the correlations with the behaviors
        mean_groom_pop = mean(rastermap(gpop,:));
        for jj = 1:size(long_eps,1)      
            try
                long_groom{i}{jj} = mean_groom_pop(long_eps(jj,1)-5*fs:long_eps(jj,2)); 
            catch
                continue
            end
        end
    
        clear allcorr mallcorr
        anym = flrthresh | fllthresh;
        for jj = 1:size(rastermap,1)
            nanidx = isnan(rastermap(jj,:));
            allcorr(jj) = corr(rastermap(jj,~nanidx)', episodes(~nanidx));
            mallcorr(jj) = corr(rastermap(jj,~nanidx)', anym(~nanidx));
            
        end
        ncorr = cat(1, ncorr, allcorr');
        gcorr = cat(1, gcorr, allcorr(gpop)');
    
        gcorrm = cat(1, gcorrm, mallcorr(gpop)');
    end

end

%% Plot histogram of neuron correlations with behavior highlight grooming population
clc
figure, gg=histogram(ncorr, 'FaceColor', [0 0 0], 'FaceAlpha', 0.3);
hold on
histogram(gcorr, 'BinWidth', gg.BinWidth, 'FaceColor', [0 0 1], 'FaceAlpha', 1)
xlabel('Correlation with Grooming')
ylabel('# Neurons')
legend({'All neurons', 'Grooming neurons', ''}, 'Location', 'Best')
gg.Parent.FontSize = 12;

% for stats comparison, remove the grooming population from rest of
% neuronal population
ncorr2 = ncorr;
for i = 1:length(gcorr)
    ncorr2(find(ncorr2==gcorr(i),1))=[];
end

h1 = isnormal(ncorr2);
h2 = isnormal(gcorr);
if ~h1 || ~h2
    disp('At least one distribution is not normal')
else
    disp('Both distributions normal')
end

disp(['Grooming neurons: ', num2str(median(gcorr)),'+/-', num2str(iqr(gcorr))])
disp(['All other neurons: ', num2str(median(ncorr2)),'+/-', num2str(iqr(ncorr2))])
[p,h, stats] = ranksum(ncorr2, gcorr);
disp(['Grooming population is different from all neurons: Wilcoxon RankSum test, p=', num2str(p)])
axis([-1 1 0 500])
% saveas(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'neuron_grooming_correlation.svg']))
% exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'neuron_grooming_correlation.png']), 'Resolution', 300)

%% Compare correlation of grooming pop to grooming vs movement

plot_data = [gcorr gcorrm];
figure, swarmchart(repmat([1 2], size(gcorrm,1), 1), plot_data, 'b', 'filled','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25, 'XJitterWidth', 0.25)
hold on
boxplot(plot_data, 'Colors', 'k', 'Symbol', '')
ylabel('Pearson correlation coefficient')
xticklabels({'Grooming', 'Movement'})
title('Grooming neuron population')
ax = gca;
ax.FontSize = 14;

h1 = isnormal(gcorr);
h2 = isnormal(gcorrm);
if ~h1 || ~h2
    disp('At least one distribution is not normal')
else
    disp('Both distributions normal')
end

[h,p, ~, stats] = ttest(gcorr, gcorrm);
disp(['Grooming: ', num2str(mean(gcorr)),'+/-', num2str(std(gcorr))])
disp(['Movement: ', num2str(mean(gcorrm)),'+/-', num2str(std(gcorrm))])
disp(['Grooming population vs all neurons: Paired t-test, p=', num2str(p)])
% saveas(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'groomingpop_groomingvsmovementcorrelation.svg']))

%% Get long grooming episodes into matrix form

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

%% Plot heatmap of long episodes
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


figure, axis off, colorbar, caxis([-1 3]), colormap(bluewhitered())
% exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'longgroom_heatmap_2p_colorbar.png']), 'Resolution', 300)
%% Plot shaded error bat for long episodes

figure, hold on
shadedErrorBar(t(1:p50), mean(long_groom_mat(:,1:p50), 'omitnan'), std(long_groom_mat(:,1:p50), 'omitnan'), 'k', 1)
line([0 0], ylim, 'Color', [0 0 0], 'LineWidth', 2)
axis([-5 30 -1 2.5])

% xticklabels([])
% yticklabels([])
% exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'longgroom_average2p.png']), 'Resolution', 300)

%% Compute and plot slope of average across long groom episodes

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

blen = 5;

clear p
for i = 1:size(long_groom_mat,1)
    
    X = long_groom_mat(i,round(blen*fs):ep_length(i));
    t = xt(X, fs);
    avg_slope(i,:) = polyfit(t, X, 1);

    % one episode has nans at end - likely due to episode occuring at the
    % end of the session
    if isnan(avg_slope(i,1))
        X = X(1:find(isnan(X),1)-1);
        t = xt(X, fs);
        avg_slope(i,:) = polyfit(t, X, 1);
    end

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
end

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
imagesc(-groom_prop), caxis([-0.5 0]); colormap(bluewhitered()), colorbar

for p = 1:58 % this is length(dorsalMaps.edgeOutline)-2 to ignore OB
    plot(dorsalMaps.edgeOutline{p}(:,2), dorsalMaps.edgeOutline{p}(:,1), 'k', 'LineWidth', 2)
end
axis equal off
set(gca, 'YDir', 'reverse');

% exportgraphics(gcf, fix_path(['Y:\nick\behavior\grooming\figures\', 'grooming_proportion2p.png']), 'Resolution', 300)



%% ----------------------------


clc, clear, close all
fileID = fopen('expt1_datalist.txt','r');
addpath('C:\Users\user\Documents\Nick\grooming\utils')

formatSpec = '%s';



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



mp_list = fix_path(mp_list);
data_list = mp_list;
current_mouse = '';
cluster_data = fix_path('Y:\nick\behavior\grooming\20241114092737_behavior_clustering.mat');
dmetric = 'cosine';

addpath(genpath('C:\Users\user\Documents\Nick\grooming\utils'))
addpath(genpath('/home/user/Documents/grooming/utils'))
try load('C:\Users\user\Documents\Nick\grooming\utils\allen_map\allenDorsalMap.mat');
catch
    load('/home/user/Documents/grooming/utils/allen_map/allenDorsalMap.mat');
end

fs = 90;
include_boris = true;


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

% trueM is the cluster results on the data averaged across all mice
dat = mean(sim_matrix,3,'omitnan');
[trueM,Q]=community_louvain(dat);

% M is the cluster results on the data for individual mice
% ari is the adjusted rand index between each mouse's clusters and averaged
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

