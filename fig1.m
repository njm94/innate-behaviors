
%% average events
addpath('C:\Users\user\Documents\Nick\ridgeModel');
addpath('C:\Users\user\Documents\Nick\ridgeModel\widefield')
addpath('C:\Users\user\Documents\Nick\ridgeModel\smallStuff') 
addpath('C:\Users\user\Documents\Nick\grooming\utils')



clc, clear
fileID = fopen('expt1_datalist.txt','r');

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = '';

fs = 90 ;
[b, a] = butter(2, 0.01/(fs/2), 'high');

thy1_idx = 1:7;
ai94_idx = 8:13;
camk_idx = 14:25;

spontaneous = [1,2,8,9,14,15,20,21];
evoked = [3:7,10:13,16:19,22:25];

states = ["Stop", "Elliptical", "Asymmetric", "Bilateral", "Unilateral"];

 

%%
N = length(data_list{1});
trial_length = get_min_trial_length_from_expt_list(data_list{1});

event_raster = zeros(N, trial_length);
num_episodes = zeros(N, 1);

tmat = zeros(numel(states), numel(states), N);

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
    [episodes, idx] = aggregate(bmat, 3);
    num_episodes(j) = size(idx,1);

    event_raster(j,:) = episodes;



%     % consolidate unilateral events
    bmat2 = bmat;
    bmat2(bmat2==3)=2;
    bmat2(bmat2==4)=3;
    bmat2(bmat==5 | bmat==6) = 4;
    bmat2(bmat==7) = 0;% 5;
    bmat2 = bmat2 + 1;

    
    
    % update transition matrix
    for ii = 1:size(idx,1)
        tmp = bmat2(idx(ii,1):idx(ii,2));
        event_idx = [1 tmp] > 1;
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



function [bmat2, idx] = aggregate(bmat, W, fs)
% aggregate grooming events by merging all events that occur within W
% seconds of each other. Then return list of start and stop indexes


if nargin < 2 || isempty(W), W = 3; end
if nargin < 3 || isempty(fs), fs = 90; end

W = round(W*fs); % convert seconds to samples

bmat2 = bmat>0;
bmat2 = movmax(bmat2, W);
bmat2 = movmin(bmat2, W);

% figure, plot(bmat>0), hold on, plot(bmat2-1)

time_on = find(diff(bmat2)==1)+1;
time_off = find(diff(bmat2)==-1)+1;

% edge cases where grooming occurs at the beginning or end
if bmat2(end), time_off = [time_off, length(bmat2)]; end
if bmat2(1), time_on = [1, time_on]; end

idx = [time_on; time_off]';
end
