
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


mp_list = {'Y:\nick\behavior\grooming\2p\ETR2_thy1\20231113143925'; ...
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231113155903'; ...
    'Y:\nick\behavior\grooming\2p\ETR3_thy1\20231115174148'};

data_list{1} = [data_list{1}; mp_list];

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
    "Elliptical Asymmetric", "Bilateral", "Unilateral"];

aggregation_sz = 3;

%%
N = length(data_list{1});

tmat = zeros(numel(states), numel(states), N);

for j = [thy1_idx, camk_idx(1), mp_idx]%1:N
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

            [y,x] = set_xy_states(last_event, current_event);

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
        [y,x] = set_xy_states(last_event, current_event);
        tmat(y,x,j) = tmat(y,x,j) + 1;

        disp([num2str(counter), ' events found'])
    end
end

%% Create conditional probability matrix from transition matrix

pmat = tmat ./ sum(tmat,2);

B = mean(pmat, 3, 'omitnan');
figure, imagesc(B)
colormap(flipud(colormap('gray'))), colorbar
xticklabels(states)
yticklabels(states)

%%
% ignore large bilateral events since they are so rare

[M,Q]=community_louvain(B([1:4,6], [1:4,6]));

%%
clc
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
pgraph = digraph(B([1:4,6], [1:4, 6]));
figure, plot(pgraph, 'MarkerSize', 15, 'LineWidth', pgraph.Edges.Weight*10, ...
    'NodeColor', cols, 'NodeFontSize', 15, ...
    'EdgeColor', 'k', 'ArrowSize', 15, 'NodeLabel', states([1:4, 6]), ...
    'Layout', 'force', 'WeightEffect', 'inverse')


%% 
function [y,x] = set_xy_states(last_state, current_state)
            switch last_state
                case 'Stationary'
                    y = 1;
                case 'Elliptical'
                    y = 2;
                case 'Right Asymmetric' 
                    y = 3;
                case 'Left Asymmetric'
                    y = 3;
                case 'Elliptical Asymmetric'
                    y = 4;
                case 'Large Bilateral'
                    y = 5;
                case 'Right' 
                    y = 6;
                case 'Left'
                    y = 6;
            end


            switch current_state
                case 'Stationary'
                    x = 1;
                case 'Elliptical'
                    x = 2;
                case 'Right Asymmetric' 
                    x = 3;
                case 'Left Asymmetric'
                    x = 3;
                case 'Elliptical Asymmetric'
                    x = 4;
                case 'Large Bilateral'
                    x = 5;
                case 'Right' 
                    x = 6;
                case'Left'
                    x = 6;
            end
end


