
%% average events
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

fs = 90 ;
% [b, a] = butter(2, 0.01/(fs/2), 'high');

thy1_idx = 1:7;
ai94_idx = 8:13;
camk_idx = 14:25;

left_idx = [6:7,13,19,25];
right_idx = [3:6,10:12,16:18,22:24];

spontaneous = [1,2,8,9,14,15,20,21];
evoked = [3:7,10:13,16:19,22:25];

states = ["Stop", "Elliptical", "Asymmetric", "Bilateral", "Unilateral"];

aggregation_sz = 3;

%%
N = length(data_list{1});
trial_length = get_min_trial_length_from_expt_list(data_list{1});

event_raster = zeros(N, trial_length);
num_episodes = zeros(N, 1);

tmat = zeros(numel(states), numel(states), N);

binned_bmat = zeros(4, 20, N);

for j = 1:N
    data_dir = data_list{1}{j};
    timestamp_file = [data_dir, filesep, getAllFiles(data_dir, '_trim.txt')];
    [mouse_root_dir, exp_date, ~] = fileparts(data_dir);
    [~, mouse_id, ~] = fileparts(mouse_root_dir);
    [expt_root_dir, ~, ~] = fileparts(mouse_root_dir);
    snippets_dir = [data_dir, filesep, 'snippets'];
%     output_dir = [data_dir, filesep, 'outputs'];

%     if isfile([output_dir, filesep, 'average_dFF.pdf'])
%         disp('Average dFF already exists. Moving on...')
%         continue
%     end
% 
%     if ~isfolder(output_dir)
%         mkdir(output_dir)
%     end


    [snippets, labels] = parse_snippets(snippets_dir);
%     num_events(j) = cellfun(@(data) size(data, 1), snippets);

    bmat = zeros(1, trial_length);
    for ii = 1:length(snippets)
        for jj = 1:size(snippets{ii},1)
            if snippets{ii}(jj,1) >= trial_length, continue; end
            bmat(snippets{ii}(jj,1):min([snippets{ii}(jj,2), trial_length])) = ii;
        end
    end
    [episodes, idx] = aggregate(bmat, aggregation_sz);
    num_episodes(j) = size(idx,1);

    event_raster(j,:) = episodes;
    all_event_idx{j} = idx;

    num_left(j) = size(snippets{matches(labels, 'left')},1);
    num_right(j) = size(snippets{matches(labels, 'right')},1); 

    num_largeleft(j) = size(snippets{matches(labels, 'largeleft')},1);
    num_largeright(j) = size(snippets{matches(labels, 'largeright')},1);

    % consolidate events to elliptical, asymmetric, bilateral, unilateral
    bmat2 = zeros(4, size(bmat,2));
    bmat2(1,:) = bmat == find(matches(labels, 'elliptical'));
    bmat2(2,:) = bmat == find(matches(labels, 'largeleft')) |  bmat == find(matches(labels, 'largeright')) ;
    bmat2(3,:) = bmat == find(matches(labels, 'largebilateral'));
    bmat2(4,:) = bmat == find(matches(labels, 'left')) |  bmat == find(matches(labels, 'right')) ;


    bmat2 = zeros(1, size(bmat,2));
    bmat2(bmat == find(matches(labels, 'elliptical'))) = 1;
    bmat2(bmat == find(matches(labels, 'largeleft')) |  bmat == find(matches(labels, 'largeright'))) = 2;
    bmat2(bmat == find(matches(labels, 'largebilateral'))) = 3;
    bmat2(bmat == find(matches(labels, 'left')) |  bmat == find(matches(labels, 'right'))) = 4;


    

%     bmat2 = bmat;
%     bmat2(bmat2==3)=2;
%     bmat2(bmat2==4)=3;
%     bmat2(bmat==5 | bmat==6) = 4;
%     bmat2(bmat==7) = 0;% 5;
    bmat2 = bmat2 + 1;

    bmat3(:,j) = bmat2(1:trial_length);


    % create bins of a minute long
    % start at end of trial and work backwards, since start is baseline
    w = 60;
    bin_edges = fliplr(trial_length:-ceil(fs*w):1);
    for ii = 1:length(bin_edges)-1
        tmp = bmat2(bin_edges(ii):bin_edges(ii+1));
        for jj = 2:5
            binned_bmat(jj-1,ii,j) = length(find(tmp==jj));
        end
    end


    
    
    % update transition matrix
    for ii = 1:size(idx,1)
        tmp = bmat2(idx(ii,1):idx(ii,2));
        event_idx = [1 tmp] > 1;
        if any(event_idx)
            event_start = find(diff(event_idx) == 1);
            for jj = 1:length(event_start)
                if jj == 1
                    tmat(1, tmp(event_start(jj)), j) = tmat(1, tmp(event_start(jj)), j)+1;
                else
                    tmat(last_event, tmp(event_start(jj)), j) = tmat(last_event, tmp(event_start(jj)), j) + 1;
                end
                last_event = tmp(event_start(jj));
            end
            tmat(last_event, 1, j) = tmat(last_event, 1, j) + 1;
        end
    end

end

%%
figure
subplot(2,2,1), hold on, 
data2plot = [mean(num_right(left_idx)), mean(num_left(left_idx))];
bar(data2plot, 0.4), 
swarmchart(ones(size(left_idx))-0.45, num_right(left_idx), 'ko', 'XJitterWidth', 0.2)
swarmchart(ones(size(left_idx))+1.45, num_left(left_idx), 'ko', 'XJitterWidth', 0.2)
ylabel('# Unilateral Strokes')
xticks([1, 2])
xticklabels({'Right', 'Left'})

[h,p] = ttest(num_right(left_idx), num_left(left_idx));
ax = gca;
ax.FontSize = 12;
title(['p=', num2str(p)])






subplot(2,2,2), hold on
data2plot = [mean(num_right(right_idx)), mean(num_left(right_idx))];
bar(data2plot, 0.4), 
swarmchart(ones(size(right_idx))-0.45, num_right(right_idx), 'ko', 'XJitterWidth', 0.2)
swarmchart(ones(size(right_idx))+1.45, num_left(right_idx), 'ko', 'XJitterWidth', 0.2)
xticks([1, 2])
xticklabels({'Right', 'Left'})

[h,p] = ttest(num_right(right_idx), num_left(right_idx));
ax = gca;
ax.FontSize = 12;
title(['p=', num2str(p)])




subplot(2,2,3), hold on
data2plot = [mean(num_largeright(left_idx)), mean(num_largeleft(left_idx))];
bar(data2plot, 0.4), 
swarmchart(ones(size(left_idx))-0.45, num_largeright(left_idx), 'ko', 'XJitterWidth', 0.2)
swarmchart(ones(size(left_idx))+1.45, num_largeleft(left_idx), 'ko', 'XJitterWidth', 0.2)
ylabel('# Asymmetric Bilateral Strokes')
xticks([1, 2])
xticklabels({'Right', 'Left'})

[h,p] = ttest(num_largeright(left_idx), num_largeleft(left_idx));
ax = gca;
ax.FontSize = 12;
title(['p=', num2str(p)])
xlabel('Left Stimulus')


subplot(2,2,4), hold on


data2plot = [mean(num_largeright(right_idx)), mean(num_largeleft(right_idx))];
bar(data2plot, 0.4), 
swarmchart(ones(size(right_idx))-0.45, num_largeright(right_idx), 'ko', 'XJitterWidth', 0.2)
swarmchart(ones(size(right_idx))+1.45, num_largeleft(right_idx), 'ko', 'XJitterWidth', 0.2)
xticks([1, 2])
xticklabels({'Right', 'Left'})

[h,p] = ttest(num_largeright(right_idx), num_largeleft(right_idx));
ax = gca;
ax.FontSize = 12;
title(['p=', num2str(p)])
xlabel('Right Stimulus')


%%

pmat = tmat ./ sum(tmat,2);
%%
B = mean(pmat(:,:,evoked), 3, 'omitnan');

%%
[M,Q]=community_louvain(B);

%%
clc
cols = zeros(length(M), 3);
col1 = [0.9290 0.6940 0.1250];
col2 = [0.3010 0.7450 0.9330];
for i = 1:length(M)
    if M(i) == 1
        cols(i,:) = col1;
    elseif M(i) == 2
        cols(i,:) = col2;
    end
end
pgraph = digraph(B);
figure, plot(pgraph, 'MarkerSize', 15, 'LineWidth', pgraph.Edges.Weight*10, ...
    'NodeColor', cols, 'NodeFontSize', 15, ...
    'EdgeColor', 'k', 'ArrowSize', 15, 'NodeLabel', states, ...
    'Layout', 'force', 'WeightEffect', 'inverse')


%%

figure, subplot(2,3,1:2)
imagesc(1-[event_raster(spontaneous,:); event_raster(evoked,:)])
colormap gray, 
axis off
y = [length(spontaneous)+0.5, length(spontaneous)+0.5, size(event_raster,1)+0.5, size(event_raster,1)+0.5];
patch([0 size(event_raster,2) size(event_raster,2) 0], y, [0.3010 0.7450 0.9330], 'FaceAlpha', 0.1, 'EdgeColor', 'none')


subplot(2,3,3)
num_spon_events = num_episodes(spontaneous);
num_evoked_events = num_episodes(evoked);
[h,p] = ttest2(num_spon_events, num_evoked_events);

size_diff = length(num_evoked_events) - length(num_spon_events);
num_spon_events = cat(1, num_spon_events, nan(size_diff, 1));
data2plot = [num_spon_events, num_evoked_events];

boxplot(data2plot, 'colors', 'kk', 'Labels',{'Spontaneous','Evoked'});
ylabel('# Grooming episodes')
hold on, swarmchart(ones(size(num_spon_events))+0.3, num_spon_events, 'ko', 'XJitterWidth', 0.25)
swarmchart(ones(size(num_evoked_events))+1.3, num_evoked_events, 'ko', 'XJitterWidth', 0.25)
ax = gca;
ax.FontSize = 12;
title(['p=', num2str(p)])


spon_probability = mean(event_raster(spontaneous, :));
spon_probability = movmean(spon_probability, round(5*fs));
evoked_probability = mean(event_raster(evoked, :));
evoked_probability = movmean(evoked_probability, round(5*fs));
t = xt(event_raster, fs, 2);


subplot(2,3,4:5)
% plot(xFit, yFit, 'k-', 'LineWidth', 1); % Plot fitted line.
hold on

plot(t, evoked_probability, 'Color', [0.3010 0.7450 0.9330],  'LineWidth', 1.5)
plot(t, spon_probability, 'k', 'LineWidth', 1.5)

drop_onset = 32:60:t(end);
vline(drop_onset, 'k:')

axis tight
ylabel('Grooming probability')
xlabel('Time (s)')
ax = gca;
ax.FontSize = 12;


subplot(2,3,6)
w = 5;
for i = 1:length(drop_onset)
    evoked_prob_centered(:,:,i) = event_raster(evoked, round(fs*(drop_onset(i)-w)): round(fs*(drop_onset(i)+ w)));
end

spon_dur = sum(event_raster(spontaneous,:),2)/fs;
evoked_dur = sum(event_raster(evoked,:),2)/fs;
spon_dur = cat(1, spon_dur, nan(size_diff, 1));
data2plot = [spon_dur, evoked_dur];

boxplot(data2plot, 'colors', 'kk', 'Labels',{'Spontaneous','Evoked'});
ylabel('Grooming duration (s)')
hold on, swarmchart(ones(size(spon_dur))+0.3, spon_dur, 'ko', 'XJitterWidth', 0.25)
swarmchart(ones(size(evoked_dur))+1.3, evoked_dur, 'ko', 'XJitterWidth', 0.25)
ax = gca;
ax.FontSize = 12;
[h,p] = ttest2(spon_dur, evoked_dur);
title(['p=', num2str(p)])


%% grooming vs days


thy1_groom = num_episodes(thy1_idx);
ai94_groom = num_episodes(ai94_idx);
camk_groom1 = num_episodes(camk_idx(1:6));
camk_groom2 = num_episodes(camk_idx(7:end));



all_groom = [thy1_groom(1:end-1), ai94_groom, camk_groom1, camk_groom2];

subplot(2,3,6),                   
plot(all_groom), hold on
boxplot(all_groom', 'colors', 'k')
xlabel('Day')
ylabel('# Grooming episodes')
ax = gca;
ax.FontSize = 12;


%% relative frequency of events for each mice spontaneous vs evoked

state_count = squeeze(sum(tmat(2:end, :, :), 2));

% compute for each mouse
thy1_relfreq = state_count(:, thy1_idx) ./ sum(state_count(:, thy1_idx));
ai94_relfreq = state_count(:, ai94_idx) ./ sum(state_count(:, ai94_idx));
camk_relfreq = state_count(:, camk_idx) ./ sum(state_count(:, camk_idx));


all_relfreq = cat(3, thy1_relfreq(:,1:6), ai94_relfreq, camk_relfreq(:,1:6), camk_relfreq(:,7:end));

all_spon_relfreq = reshape(all_relfreq(:,1:2,:), 4, []);
all_evoked_relfreq = reshape(all_relfreq(:,3:6,:), 4, []);


% % compute on compiled data from all mice
all_relfreq = state_count ./ sum(state_count);
all_spon_relfreq = sum(state_count(:, spontaneous), 2)./ sum(state_count(:, spontaneous), [2 1]);
all_evoked_relfreq = sum(state_count(:, evoked), 2) ./ sum(state_count(:, evoked), [2 1]);



figure, 
b=bar([mean(all_spon_relfreq,2), mean(all_evoked_relfreq,2)]);
b(1).FaceColor = 'k';
b(2).FaceColor =  [0.3010 0.7450 0.9330];

hold on
% errorbar((1:4)-0.15, mean(all_spon_relfreq,2), std(all_spon_relfreq,[],2)./sqrt(size(all_spon_relfreq,2)), 'k.')
% errorbar((1:4)+0.15, mean(all_evoked_relfreq,2), std(all_evoked_relfreq,[],2)./sqrt(size(all_evoked_relfreq,2)), 'k.')
xticklabels(states(2:end))
ylabel('Relative Frequency')
legend(["Spontaneous", "Evoked"], 'Location', 'Best', 'Box', 'off')
ax = gca;
ax.FontSize = 12;

% swarmchart((1:4)-0.35, all_spon_relfreq, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [1 1 1], 'XJitter', 'density');
% swarmchart((1:4)+0.35, all_evoked_relfreq, 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerEdgeColor', [0.3010 0.7450 0.9330])


%%

summary_binned_bmat = 100*(sum(binned_bmat(:,:,evoked),3)./sum(binned_bmat(:,:,evoked), [3 1]))'; 

sliding_bmat_for_area = squeeze(mean(sliding_bmat(:, evoked, :), 2, 'omitnan'));
figure, 
% area(summary_binned_bmat(:,[1 2 4 3]))
area(t, sliding_bmat_for_area)
% legend(circshift(["Elliptical", "Asymmetric", "Bilateral", "Unilateral"], 1))
xlabel('Time (mins)')
ylabel('Proportion of behavior')
axis([1 t(end) 0 1])
ax = gca;
ax.FontSize = 12;

%%
sliding_bmat = zeros(trial_length, N, 4);

w = 300;

for i = 2:5
    sliding_bmat(:,:,i-1) = movsum(bmat3==i,  round(w*fs), 1) ./ movsum(bmat3~=1, round(w*fs), 1);
end
avg_sliding_bmat_thy1 = squeeze(mean(sliding_bmat(:, thy1_idx(3:end), :), 2, 'omitnan'));
avg_sliding_bmat_thy1 = 100 * (avg_sliding_bmat_thy1 - avg_sliding_bmat_thy1(1,:));
avg_sliding_bmat_ai94 = squeeze(mean(sliding_bmat(:, ai94_idx(3:end), :), 2, 'omitnan'));
avg_sliding_bmat_ai94 = 100 * (avg_sliding_bmat_ai94 - avg_sliding_bmat_ai94(1,:));
avg_sliding_bmat_camk1 = squeeze(mean(sliding_bmat(:, camk_idx(3:6), :), 2, 'omitnan'));
avg_sliding_bmat_camk1 = 100 * (avg_sliding_bmat_camk1 - avg_sliding_bmat_camk1(1,:));
avg_sliding_bmat_camk2 = squeeze(mean(sliding_bmat(:, camk_idx(9:end), :), 2, 'omitnan'));
avg_sliding_bmat_camk2 = 100 * (avg_sliding_bmat_camk2 - avg_sliding_bmat_camk2(1,:));


avg_sliding_bmat = mean(cat(3, avg_sliding_bmat_thy1, avg_sliding_bmat_ai94, avg_sliding_bmat_camk1, avg_sliding_bmat_camk2), 3);
sem_sliding_bmat = std(cat(3, avg_sliding_bmat_thy1, avg_sliding_bmat_ai94, avg_sliding_bmat_camk1, avg_sliding_bmat_camk2), [], 3) ./sqrt(4);

% avg_sliding_bmat = squeeze(mean(sliding_bmat(:, evoked, :), 2, 'omitnan'));

figure, hold on

% plot(t, avg_sliding_bmat)
cols = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"];
for i = 1:4
    s = shadedErrorBar(t, avg_sliding_bmat(:,i), sem_sliding_bmat(:,i), [], 1);
    set(s.edge,'Color', cols(i))
    s.mainLine.Color = cols(i);
    s.mainLine.LineWidth = 2;
    s.patch.FaceColor = cols(i);
end

% plot(t, 100*(avg_sliding_bmat_thy1(:,4) - avg_sliding_bmat_thy1(1,4)))
% plot(t, 100*(avg_sliding_bmat_ai94(:,4) - avg_sliding_bmat_ai94(1,4)))
% plot(t, 100*(avg_sliding_bmat_camk1(:,4) - avg_sliding_bmat_camk1(1,4)))
% plot(t, 100*(avg_sliding_bmat_camk2(:,4) - avg_sliding_bmat_camk2(1,4)))

ylabel('\Delta % from T_0')
xlabel('Time (s)')
title([num2str(w/60), ' min sliding window'])
% legend(["thy1", "ai94", "hyl3", "ibl2"])
legend(["", "", "", "Elliptical", "", "", "", "Asymmetric", "", "", "", "Bilateral", "", "", "", "Unilateral"], 'Location', 'Best')
hline(0, 'k:')

%%
sx = [];
sy = [];
figure, hold on
for i = 1:length(evoked)
    idx = all_event_idx{evoked(i)};
    sx = cat(1, sx, idx(:,1));
    sy = cat(1, sy, diff(idx, 1, 2));
    
end

sx = sx/fs;
sy = sy/fs;

scatter(sx, sy, 'filled')
lsline

%%

binsx = 1:60:max(sx);
for i = 1:length(binsx)-1
    tmp = find(sx>binsx(i) & sx <= binsx(i+1));
    avg_dur(i) = mean(sy(tmp));
    num_episodes(i) = numel(tmp);
end
figure, 
subplot(1,2,1),
plot(avg_dur)
ylabel('Average Duration of Grooming Episode')
xlabel('Time (mins)')
subplot(1,2,2)
plot(num_episodes)
ylabel('Number of grooming episodes')
xlabel('Time (mins)')



%%





%%



function min_trial_length = get_min_trial_length_from_expt_list(expt_list)
    N = length(expt_list);
    trial_length = zeros(1,N);
    
    for j = 1:N
        data_dir = expt_list{j};
        timestamp_file = [data_dir, filesep, getAllFiles(data_dir, '_trim.txt')];

        % take 5 frames off the end to make sure we take out the LED OFF
        timestamps = readmatrix(timestamp_file);
        trial_length(j) = timestamps(end-5,1);
    end

    min_trial_length = min(trial_length);
end




