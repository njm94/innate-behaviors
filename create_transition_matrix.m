
%% 
addpath('C:\Users\user\Documents\Nick\ridgeModel');
addpath('C:\Users\user\Documents\Nick\ridgeModel\widefield')
addpath('C:\Users\user\Documents\Nick\ridgeModel\smallStuff') 
addpath('C:\Users\user\Documents\Nick\grooming\utils')
addpath('C:\Users\user\Documents\Nick\grooming\utils\layoutCode')


clc, clear
fileID = fopen('expt1_datalist.txt','r');

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = ''; 


% mp_list = {'Y:\nick\behavior\grooming\2p\ETR2_thy1\20231113143925'; ...
%     'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231113155903'; ...
%     'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231115174148'};
mp_list = {'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240729'; ...
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240731';
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240802';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240731';
    'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240729';
    'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240802'};


% data_list{1} = mp_list;

fs = 90 ;

thy1_idx = 1:7;
ai94_idx = 8:13;
camk_idx = 14:25;
mp_idx = 26:28;

left_idx = [6:7,13,19,25];
right_idx = [3:6,10:12,16:18,22:24];

spontaneous = [1,2,8,9,14,15,20,21];
evoked = [3:7,10:13,16:19,22:25];

states = ["Stationary", "Elliptical", "Asymmetric", ...
    "Bilateral", "Elliptical Asymmetric", "Unilateral"];

aggregation_sz = 3;

% data_list{1} = [data_list{1}(evoked); mp_list];
data_list{1} = [data_list{1}; mp_list];
% data_list{1} = mp_list;

%%
N = length(data_list{1});

tmat = zeros(numel(states), numel(states), N);

for j = 1:N
    data_dir = data_list{1}{j};
    timestamp_file = [data_dir, filesep, getAllFiles(data_dir, '_trim.txt')];
    [mouse_root_dir, exp_date, ~] = fileparts(data_dir);
    [~, mouse_id, ~] = fileparts(mouse_root_dir);
    [expt_root_dir, ~, ~] = fileparts(mouse_root_dir);
    snippets_dir = [data_dir, filesep, 'snippets'];


    if ~isempty(getAllFiles(data_dir, '.tsv'))
        boris_file = [data_dir, filesep, getAllFiles(data_dir, '.tsv')];
        [events, snippets, b_table] = read_boris(boris_file);

        % consolidate lick events
        lick_idx = contains(events.Properties.VariableNames, 'Lick');
        lick_events = events(:,lick_idx);
        lick_events = any(table2array(lick_events),2);
        
        % remove lick and point events from event matrix
        idx = contains(events.Properties.VariableNames, 'Lick') | ...
            contains(events.Properties.VariableNames, 'Drop') | ...
            contains(events.Properties.VariableNames, 'Video') | ...
            contains(events.Properties.VariableNames, 'Flail') ;
        stroke_events = removevars(events, idx);
        snippets(idx) = [];
        labels = stroke_events.Properties.VariableNames;

        bmat = any(table2array(stroke_events),2);
    else 
        continue
        % the snippets are a soon-to-be deprecated analysis method. Remove 
        % this when fully transitioned to the BORIS labels
        [snippets, labels] = parse_snippets(snippets_dir);
    end


    [episodes, idx] = aggregate(bmat, aggregation_sz);
    num_episodes(j) = size(idx,1);

    % concatenate all snippets into one array
    tmp_snippets = catcell(1, snippets);
    
    % create a labels matrix of same size as snippets to keep track of
    % which behavior is happening
    tmp_labels = cell(1, length(labels));
    for jj = 1:length(labels)
        tmp_labels{jj} = repmat(labels(jj), size(snippets{jj},1), 1);
    end
    tmp_labels = catcell(1, tmp_labels);

    for ii = 1:size(idx,1)
        disp(['EPISODE ', num2str(ii)])
        % define last event to be STATIONARY event
        last_event = 'Stationary';
        last_event_idx = idx(ii,1);

        counter = 0;

        while last_event_idx < idx(ii, 2) || last_event_idx == idx(ii,1)
%             disp(size(tmp_snippets))
            [~, tmp_idx] = min(abs(tmp_snippets(:,1)-last_event_idx));       
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
        current_event = 'Stationary';
        [y,x] = set_xy_states(last_event, current_event, states);
        tmat(y,x,j) = tmat(y,x,j) + 1;

        disp([num2str(counter), ' events found'])
    end
end

%% Create conditional probability matrix from transition matrix

figure,

subplot(1,2,1), imagesc(mean(tmat,3))
colormap(flipud(colormap('gray'))), 
c=colorbar;
c.Label.String = 'Count';
xticks(1:length(states))
yticks(1:length(states))
xticklabels(states)
yticklabels(states)
title('Transitions', 'FontSize', 16)
ylabel('From', 'FontSize', 14)
xlabel('To', 'FontSize', 14)


pmat = tmat ./ sum(tmat,2);

B = median(pmat, 3, 'omitnan');
subplot(1,2,2), imagesc(B)
colormap(flipud(colormap('gray'))), 
c=colorbar;
c.Label.String = 'Probability';
title('Transition probability', 'FontSize', 16)
xticks(1:length(states))
yticks(1:length(states))
xticklabels(states)
yticklabels(states)
ylabel('From', 'FontSize', 14)
xlabel('To', 'FontSize', 14)
%% For Louvain Clustering, eliminate stationary from group by introducing recursive connection instead
new_tmat = tmat(2:end, 2:end, :);
for i = 1:size(tmat,3)
    for j = 1:size(new_tmat,2)
        num_stationary_to_event_j = tmat(1,j+1,i);
        num_event_j_to_stationary = tmat(j+1, 1, i);
        new_tmat(j,j,i) = new_tmat(j,j,i) + num_stationary_to_event_j + num_event_j_to_stationary;

    end
end



%%


p1 = nanmean(pmat(:,:,thy1_idx), 3);
p2 = nanmean(pmat(:,:,camk_idx(1:5)), 3);
p3 = nanmean(pmat(:,:,mp_idx), 3);



figure
subplot(1,3,1)
imagesc(p1), caxis([0 0.6])
colormap(flipud(colormap('gray'))), colorbar
xticks(1:6)
xticklabels(states)
yticklabels(states)

subplot(1,3,2)
imagesc(p2), caxis([0 0.6])
colormap(flipud(colormap('gray'))), colorbar
xticks(1:6)
xticklabels(states)
yticklabels(states)

subplot(1,3,3)
imagesc(p3), caxis([0 0.6])
colormap(flipud(colormap('gray'))), colorbar
xticks(1:6)
xticklabels(states)
yticklabels(states)
%%
% ignore large bilateral events since they are so rare

test_gamma = 0.1:0.1:3;
for i = 1:length(test_gamma)
    [M,Q(i)]=community_louvain(B, test_gamma(i));
    num_uniq(i) = numel(unique(M));
end

figure, plot(test_gamma, Q);
hold on
plot(test_gamma, num_uniq)
%%
% dat = B(2:6, 2:6, 1);
dat = B;

[M,Q]=community_louvain(dat, 0.9);
% [M,Q]=community_louvain(btest);

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
pgraph = digraph(dat);
% pgraph = digraph(B(2:6, 2:6));
% pgraph = digraph(B);
% pgraph = digraph(btest);
figure, plot(pgraph, 'MarkerSize', 15, 'LineWidth', pgraph.Edges.Weight*10, ...
    'NodeColor', cols, 'NodeFontSize', 15, ...
    'EdgeColor', 'k', 'ArrowSize', 15, ...
    'NodeLabel',states,...states(2:6), ... states,
    'Layout', 'force3', 'WeightEffect', 'inverse')
axis off

%%
dat = pmat;
num_behaviors = length(unique(pgraph.Edges.EndNodes));
weights = zeros(num_behaviors, num_behaviors-1, size(dat,3));
for i = 1:num_behaviors
    for j = 1:size(pgraph.Edges.EndNodes, 1)
    end
end


%%
dat = B(2:end, 2:end);
clear tmp
for i = 1:size(dat,1)
    test = dat(i,:);
    test(i) = [];

    tmp(:,i) = test;
end



%% 
function [y,x] = set_xy_states(last_state, current_state, state_order)
        stat_idx = find(strcmpi(state_order, 'Stationary'));
        ellip_idx = find(strcmpi(state_order, 'Elliptical'));
        assym_idx = find(strcmpi(state_order, 'Asymmetric'));
        bilat_idx = find(strcmpi(state_order, 'Bilateral'));
        ellip_assym_idx = find(strcmpi(state_order, 'Elliptical Asymmetric'));
        uni_idx = find(strcmpi(state_order, 'Unilateral'));

            switch last_state
                case 'Stationary'
                    y = stat_idx;
                case 'Elliptical'
                    y = ellip_idx;
                case 'Right Asymmetric' 
                    y = assym_idx;
                case 'Left Asymmetric'
                    y = assym_idx;
                case 'Elliptical Asymmetric'
                    y = ellip_assym_idx;
                case 'Large Bilateral'
                    y = bilat_idx;
                case 'Right' 
                    y = uni_idx;
                case 'Left'
                    y = uni_idx;
            end


            switch current_state
                case 'Stationary'
                    x = stat_idx;
                case 'Elliptical'
                    x = ellip_idx;
                case 'Right Asymmetric' 
                    x = assym_idx;
                case 'Left Asymmetric'
                    x = assym_idx;
                case 'Elliptical Asymmetric'
                    x = ellip_assym_idx;
                case 'Large Bilateral'
                    x = bilat_idx;
                case 'Right' 
                    x = uni_idx;
                case'Left'
                    x = uni_idx;
            end
end


