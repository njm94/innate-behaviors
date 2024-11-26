
%% 
addpath('C:\Users\user\Documents\Nick\ridgeModel');
addpath('C:\Users\user\Documents\Nick\ridgeModel\widefield')
addpath('C:\Users\user\Documents\Nick\ridgeModel\smallStuff') 
addpath('C:\Users\user\Documents\Nick\grooming\utils')
addpath('C:\Users\user\Documents\Nick\grooming\utils\layoutCode')


clc, clear, close all
fileID = fopen('expt1_datalist.txt','r');

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = ''; 


% mp_list = {'Y:\nick\behavior\grooming\2p\ETR2_thy1\20231113143925'; ...
%     'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231113155903'; ...
%     'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231115174148'};
mp_list = {'Y:\nick\behavior\grooming\2p\ETR2_thy1\20231113143925'; ...
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231113155903'; ...
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231115174148'; ...
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240729'; ...
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240731';
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240802';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240729';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240731';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240802';
    'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240729';
    'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240802'};

cluster_data = fix_path('Y:\nick\behavior\grooming\20241114092737_behavior_clustering.mat');
% data_list{1} = mp_list;

fs = 90 ;

thy1_idx = 1:7;
ai94_idx = 8:13;
hyl3_idx = 14:19;
ibl2_idx = 20:25;
camk_idx = 14:25;
etr2_idx = 26;
etr3_idx = 27:28;
ecl3_idx = 29:31;
idr3_idx = 32:34;
rr3_idx = 35:36;

left_idx = [6:7,13,19,25];
right_idx = [3:6,10:12,16:18,22:24];

spontaneous = [1,2,8,9,14,15,20,21];
evoked = [3:7,10:13,16:19,22:25];

% states = ["Stationary", "Elliptical", "Asymmetric", ...
%     "Elliptical Right", "Elliptical Left", "Unilateral"];

% states = ["Stationary", "Elliptical", "Right Asymmetric", "Left Asymmetric" ...
%     "Elliptical Right", "Elliptical Left", "Unilateral"];

states = ["Start", "Elliptical", "Right Asymmetric", "Left Asymmetric" ...
    "Elliptical Right", "Elliptical Left", "Right", "Left", "Stop"];
% states = states(randperm(length(states)));
states = ["Start", "Right", "Left","Elliptical", "Right Asymmetric", "Left Asymmetric" ...
"Elliptical Right", "Elliptical Left", "Stop"];

% states = ["Start", "Unilateral","Elliptical", "Right Asymmetric", "Left Asymmetric" ...
% "Elliptical Right", "Elliptical Left", "Stop"];

aggregation_sz = 3;

% data_list{1} = [data_list{1}(evoked); mp_list];
data_list{1} = [data_list{1}; mp_list];
% data_list{1} = mp_list;

%%
N = length(data_list{1});

tmat = zeros(numel(states), numel(states), N);
episode_durations = cell(1, N);

eth_example = 22;
eth_counter = 1;
for j = 1:N
    data_dir = fix_path(data_list{1}{j});
    timestamp_file = [data_dir, filesep, getAllFiles(data_dir, '_trim.txt')];
    [mouse_root_dir, exp_date, ~] = fileparts(data_dir);
    [~, mouse_id, ~] = fileparts(mouse_root_dir);
    [expt_root_dir, ~, ~] = fileparts(mouse_root_dir);
    snippets_dir = [data_dir, filesep, 'snippets'];

    
    if ~isempty(getAllFiles(data_dir, '.tsv'))
        boris_file = [data_dir, filesep, getAllFiles(data_dir, '.tsv')];
        [events, snippets, b_table, video_end, cluster_labels] = get_labels_from_clustering_results(cluster_data, boris_file);
        if isempty(events) % no grooming behaviors
            continue
        end
        labels = events.Properties.VariableNames;
        bmat = any(table2array(events),2);
 
    else 
        continue
        % the snippets are a soon-to-be deprecated analysis method. Remove 
        % this when fully transitioned to the BORIS labels
        [snippets, labels] = parse_snippets(snippets_dir);
    end


    [episodes, idx] = aggregate(bmat, aggregation_sz);
    episode_durations{j} = diff(idx, 1, 2)/90;
    num_episodes(j) = size(idx,1);

    % concatenate all snippets into one array
    tmp_snippets = catcell(1, snippets);
    
    % create a labels matrix of same size as snippets to keep track of
    % which behavior is happening
%     tmp_labels = cell(1, length(labels));
%     for jj = 1:length(labels)
%         tmp_labels{jj} = repmat(labels(jj), size(snippets{jj},1), 1);
%     end
%     tmp_labels = catcell(1, tmp_labels);
    tmp_labels = cluster_labels;

    for ii = 1:size(idx,1)
        disp(['EPISODE ', num2str(ii)])
        % define last event to be STATIONARY event
        last_event = 'Start';
        last_event_idx = idx(ii,1);

        counter = 0;

        while last_event_idx < idx(ii, 2) || last_event_idx == idx(ii,1)
%             disp(size(tmp_snippets))
            [~, tmp_idx] = min(abs(tmp_snippets(:,1)-last_event_idx));       
            if isempty(tmp_idx)
                break
            end
            current_event = tmp_labels{tmp_idx};

            [y,x] = set_xy_states(last_event, current_event, states);

            tmat(y,x,j) = tmat(y,x,j) + 1;
            disp([last_event, '   ', current_event]);

            last_event = current_event;
            last_event_idx = tmp_snippets(tmp_idx, 2);

            tmp_snippets(tmp_idx,:) = [];
            tmp_labels(tmp_idx) = [];

            counter = counter + 1;
%             disp(size(tmp_snippets));
        end
        current_event = 'Stop';
        [y,x] = set_xy_states(last_event, current_event, states);
        tmat(y,x,j) = tmat(y,x,j) + 1;

        disp([num2str(counter), ' events found'])
        if counter > 50
            if eth_counter == eth_example
                data2plot = events(idx(ii,1):idx(ii,2),:);
                plot_ethogram(data2plot, states, fs)
            end
            eth_counter = eth_counter+1;
        end
    end
end

%% save the ethogram example

ax = gcf;
saveas(ax, fix_path(['Y:\nick\behavior\grooming\figures\','ethogram', '.svg']))


%%

episode_dur = catcell(1, episode_durations);

figure, 
histogram(episode_dur, 50)
set(gca, 'YScale', 'log')
ax = gcf;
saveas(ax, fix_path(['Y:\nick\behavior\grooming\figures\','episode_duration', '.svg']))


%% Create conditional probability matrix from transition matrix

figure('Position', [161 368 1261 495])
% figure('Position', [821 363 606 488])


subplot(1,2,1), 
imagesc(mean(tmat,3, 'omitnan'))
colormap(flipud(colormap('gray'))), 
c=colorbar;
c.Label.String = 'Count';
c.FontSize = 12;
xticks(1:length(states))
yticks(1:length(states))
xticklabels([])
yticklabels([])
xticklabels(states)
yticklabels(states)
title('Transitions A->B', 'FontSize', 16)
ylabel('State A', 'FontSize', 14)
xlabel('State B', 'FontSize', 14)
% axis off
box on


pmat = tmat ./ sum(tmat,2);
pmat(isnan(pmat)) = 0;
B = mean(pmat, 3, 'omitnan');
subplot(1,2,2), imagesc(B)
colormap(flipud(colormap('gray'))), 
c=colorbar;
c.Label.String = 'Probability';
c.FontSize = 12;
title('Transition probability P(B|A)', 'FontSize', 16)
xticks(1:length(states))
yticks(1:length(states))
xticklabels(states)
yticklabels(states)
ylabel('State A', 'FontSize', 14)
xlabel('State B', 'FontSize', 14)

ax = gcf;
% saveas(ax, fix_path(['Y:\nick\behavior\grooming\figures\','TPmatrices', '.pdf']))
exportgraphics(ax, fix_path(['Y:\nick\behavior\grooming\figures\','TPmatrices', '.pdf']), 'ContentType', 'vector')
% % axis off
% exportgraphics(ax, fix_path(['Y:\nick\behavior\grooming\figures\','TPmatrices', '.png']), 'Resolution', 300)


%%
% dat = B(2:6, 2:6, 1);
% dat = B;
dat = mean(pmat,3, 'omitnan');

[M,Q]=community_louvain(dat);

cols = zeros(length(M), 3);
col1 = [0 0.4470 0.7410];
col2 = [0.8500 0.3250 0.0980];
col3 = [0.9290 0.6940 0.1250];
col4 = [0.4940 0.1840 0.5560];
for i = 1:length(M)
    if M(i) == 1
        cols(i,:) = col1;
    elseif M(i) == 2
        cols(i,:) = col2;
    elseif M(i) == 3
        cols(i,:) = col3;
    else
        cols(i,:) = col4;
    end
end
pgraph = digraph(dat, 'omitselfloops');

figure('Position', [843 70 649 826]), 
plot(pgraph, 'MarkerSize', 20, 'LineWidth', pgraph.Edges.Weight*15, ...
    'NodeColor', cols, ...'NodeFontSize', 15, ...
    'EdgeColor', 'k', 'EdgeAlpha', 1, 'ArrowSize', 15, ...
    'NodeLabel','', ...
    'Layout', 'layered', 'Sources', 1, 'Sinks', length(states));%, 'WeightEffect', 'inverse')
% 
ax = gca;
% exportgraphics(ax, fix_path(['Y:\nick\behavior\grooming\figures\','network', '.png']), 'Resolution', 300)
saveas(ax, fix_path(['Y:\nick\behavior\grooming\figures\','network', '.svg']))


trueM = M;
%%


p1 = nanmean(pmat(:,:,thy1_idx), 3);
p2 = nanmean(pmat(:,:,ai94_idx), 3);
p3 = nanmean(pmat(:,:,hyl3_idx), 3);
p4 = nanmean(pmat(:,:,ibl2_idx), 3);
p5 = nanmean(pmat(:,:,etr2_idx), 3);
p6 = nanmean(pmat(:,:,etr3_idx), 3);
p7 = nanmean(pmat(:,:,ecl3_idx), 3);
p8 = nanmean(pmat(:,:,idr3_idx), 3);



figure
idat = {'p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8'};
idat = {'p4', 'p5', 'p6'};
for i = 1:length(idat)
    tmp =eval(idat{i});
    subplot(2,length(idat),i)
    imagesc(tmp), caxis([0 0.6])
    colormap(flipud(colormap('gray'))), colorbar
    xticks(1:size(tmp,2))
    yticks(1:size(tmp,1))
    xticklabels(states)
    yticklabels(states)


    [M,Q]=community_louvain(tmp);
    iM(:,i) = M;

    cols = zeros(length(M), 3);
    col1 = [0 0.4470 0.7410];
    col2 = [0.8500 0.3250 0.0980];
    col3 = [0.9290 0.6940 0.1250];
    col4 = [0.4940 0.1840 0.5560];
    for ii = 1:length(M)
        if M(ii) == 1
            cols(ii,:) = col1;
        elseif M(ii) == 2
            cols(ii,:) = col2;
        elseif M(ii) == 3
            cols(ii,:) = col3;
        else
            cols(ii,:) = col4;
        end
    end
    pgraph = digraph(tmp, 'omitselfloops');
    
    % figure('Position', [843 70 649 826]), 
    subplot(2,length(idat),i+length(idat))
    plot(pgraph, 'MarkerSize', 10, ...'LineWidth', pgraph.Edges.Weight*15, ...
        'NodeColor', cols, 'NodeFontSize', 10, ...
        'EdgeColor', 'k', 'EdgeAlpha', 1, 'ArrowSize', 10, ...
        'NodeLabel',states, ...
        'Layout', 'layered', 'Sources', 1, 'Sinks', length(states));%, 'WeightEffect', 'inverse')

    ri(i) = rand_index(trueM, iM(:,i));

end

%%

figure, boxplot(ri', 'color', 'k'), hold on
swarmchart(ones(size(ri))+0.2, ri, 'k', 'XJitterWidth', 0.2)
axis([0.5 1.75 -1 1])


%% 
function [y,x] = set_xy_states(last_state, current_state, state_order)
        start_idx = find(strcmpi(state_order, 'Start'));
        ellip_idx = find(strcmpi(state_order, 'Elliptical'));
        % assym_idx = find(strcmpi(state_order, 'Asymmetric'));
        right_assym_idx = find(strcmpi(state_order, 'Right Asymmetric'));
        left_assym_idx = find(strcmpi(state_order, 'Left Asymmetric'));
        bilat_idx = find(strcmpi(state_order, 'Bilateral'));
        ellip_right_idx = find(strcmpi(state_order, 'Elliptical Right'));
        ellip_left_idx = find(strcmpi(state_order, 'Elliptical Left'));
        % ellip_assym_idx = find(strcmpi(state_order, 'Elliptical Asymmetric'));
        uni_idx = find(strcmpi(state_order, 'Unilateral'));
        right_idx = find(strcmpi(state_order, 'Right'));
        left_idx = find(strcmpi(state_order, 'Left'));
        stop_idx =  find(strcmpi(state_order, 'Stop'));

            switch last_state
                case 'Start'
                    y = start_idx;
                case 'Elliptical'
                    y = ellip_idx;
                case 'Right Asymmetric' 
                    % y = assym_idx;
                    y = right_assym_idx;
                case 'Left Asymmetric'
                    % y = assym_idx;
                    y = left_assym_idx;
                case 'Elliptical Right'
                    % y = ellip_assym_idx;
                    y = ellip_right_idx;
                case 'Elliptical Left'
                    % y = ellip_assym_idx;
                    y = ellip_left_idx;
                case 'Large Bilateral'
                    y = bilat_idx;
                case 'Right' 
                    % y = uni_idx;
                    y = right_idx;
                case 'Left'
                    % y = uni_idx;
                    y = left_idx;
                case 'Stop'
                    y = stop_idx;
            end


            switch current_state
                case 'Start'
                    x = start_idx;
                case 'Elliptical'
                    x = ellip_idx;
                case 'Right Asymmetric' 
                    % x = assym_idx;
                    x = right_assym_idx;
                case 'Left Asymmetric'
                    % x = assym_idx;
                    x = left_assym_idx;
                case 'Elliptical Left'
                    % x = ellip_assym_idx;
                    x = ellip_left_idx;
                case 'Elliptical Right'
                    % x = ellip_assym_idx;
                    x = ellip_right_idx;
                case 'Large Bilateral'
                    x = bilat_idx;
                case 'Right' 
                    % x = uni_idx;
                    x = right_idx;
                case 'Left'
                    % x = uni_idx;
                    x = left_idx;
                case 'Stop'
                    x = stop_idx;
            end
end



function crossings = countEdgeCrossings(G, nodeOrder)
    edges = G.Edges.EndNodes;
    crossings = 0;
    for i = 1:height(edges)
        for j = i+1:height(edges)
            % Get node indices for the two edges
            edge1 = edges(i, :);
            edge2 = edges(j, :);
            % Check if the edges cross based on node order
            if (nodeOrder(edge1(1)) < nodeOrder(edge2(1)) && nodeOrder(edge1(2)) > nodeOrder(edge2(2))) || ...
               (nodeOrder(edge1(1)) > nodeOrder(edge2(1)) && nodeOrder(edge1(2)) < nodeOrder(edge2(2)))
                crossings = crossings + 1;
            end
        end
    end
end
