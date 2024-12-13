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
    % 'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231113155903';
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231115174148';};
mp_list = fix_path(mp_list);
current_mouse = '';
cluster_data = fix_path('Y:\nick\behavior\grooming\20241114092737_behavior_clustering.mat');

fs = 90;
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
    flrthresh = flrv>mean(flrv) + std(flrv);
    fllthresh = fllv>mean(fllv) + std(fllv);
    
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
            % 52 115;
            % 178 216;
            ];
    elseif i == 2
        pops = [1 32;
            70 99];
        zoom_x = [45725 88684];
    elseif i == 3
        pops = [300 380;
            385 419;
            ];
            % 1 30];
    elseif i == 4
        pops = [144 199;
            28 71];
    % elseif i == 5
    %     pops = [210 240;
    %         1 105];
    elseif i == 5
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
    for jj = 1:size(pops,1)
        line([0 0], pops(jj,:), 'Color', cols(jj,:), 'LineWidth', 5)
        plot(30*mean(rastermap(pops(jj,1):pops(jj,2),:))- 100*jj, 'LineWidth', 1, 'Color', cols(jj,:))
    end
    line([1 fs*60], [-500 -500], 'Color', [0 0 0], 'LineWidth', 2)

    % vline(zoom_x, 'k-')
    % exportgraphics(h2, fix_path(['Y:\nick\behavior\grooming\figures\thy1_mouse_ethogram_raster_zoom.png']), 'Resolution', 300)
    
    % axis([zoom_x, ylim])
    % line([zoom_x(1)+1 zoom_x(1)+1+fs*20], [-500 -500], 'Color', [0 0 0], 'LineWidth', 2)
    % axis([64253    75989 ylim])
    % line([64253 64253 + 90*10], [-50 -50], 'Color', [0 0 0], 'LineWidth', 2)
    % exportgraphics(h2, fix_path(['Y:\nick\behavior\grooming\figures\thy1_ethogram_raster_zoom.png']), 'Resolution', 300)
    
    uiopen([mouse_root, filesep, 'dalsa', filesep, 'atlas_aligned.fig'],1), hold on
    scatter(x3, y3,25, '.', 'MarkerEdgeColor', [0 0 0])
    for jj = 1:size(pops,1)
        idx = isort+1;
        scatter(x3(idx(pops(jj,1):pops(jj,2))), y3(idx(pops(jj,1):pops(jj,2))), 36, 'MarkerFaceColor', cols(jj,:), 'MarkerEdgeColor', 'k')
    end

end













%%

% figure
groom_prop = zeros(size(dorsalMaps.dorsalMapScaled));
move_prop = zeros(size(dorsalMaps.dorsalMapScaled));
none_prop = zeros(size(dorsalMaps.dorsalMapScaled));
mp = round(size(dorsalMaps.dorsalMapScaled,2)./2);
for j = 1:length(region_locs)
    % If the value is negative, it is on the right side
    addmask = dorsalMaps.dorsalMapScaled == abs(areanames.(region_locs{j}));
    if areanames.(region_locs{j}) < 0
        addmask(:, 1:mp) = 0;
    else
        addmask(:,mp:end) = 0;
    end
    groom_prop = groom_prop + addmask.* (region_neurons(j,1) ./ sum(region_neurons(j,1:3)));
    move_prop = move_prop + addmask.* (region_neurons(j,2) ./ sum(region_neurons(j,1:3)));
    none_prop = none_prop + addmask.* (region_neurons(j,3) ./ sum(region_neurons(j,1:3)));
end
%%
parent_region = cellfun(@(x) x(1:end-3), region_locs, 'UniformOutput', false);
unique_region = unique(parent_region);
clear consolidated_region_neurons
for i = 1:length(unique_region)
    idx = contains(region_locs, unique_region{i});
    consolidated_region_neurons(i,:) = sum(region_neurons(idx,:), 1);
end

%%

cols = [0 0 1; 
    1 0 0;
    0 0 0];
figure, 
h = bar((consolidated_region_neurons ./ sum(consolidated_region_neurons,2)), 'stacked');
% Apply colors
for i = 1:length(h)
    h(i).FaceColor = 'flat';       % Enable flat coloring
    h(i).CData = repmat(cols(i, :), size(consolidated_region_neurons, 1), 1); % Assign colors
end
legend({'Groom + Move', 'Move only', 'Neither'})
ylabel('Proportion')
xticklabels(unique_region)
% xlabel('Episode duration (s)')
% ax = gcf;


%%
figure, 
for i = 1:3
    subplot(1,3,i), hold on
    switch i
        case 1, imagesc(-groom_prop), caxis([-0.5 0]); colormap(bluewhitered()),  freezeColors
        case 2, imagesc(move_prop), colormap(bluewhitered()), caxis([0 0.5]); freezeColors
        case 3, imagesc(none_prop), colormap(flipud(colormap(gray))), caxis([0 1]), %h3=colorbar;
    end
    for p = 1:58%length(dorsalMaps.edgeOutline)
        % if p == 60, continue; end
        plot(dorsalMaps.edgeOutline{p}(:,2), dorsalMaps.edgeOutline{p}(:,1), 'k')
    end
            axis equal off
        set(gca, 'YDir', 'reverse');
end

%%
