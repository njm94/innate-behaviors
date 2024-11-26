%% quantify behaviors

%%
clear, clc
fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
data_list = textscan(fileID, formatSpec);

fs = 90;
% [b, a] = butter(2, 0.01/(fs/2), 'high');

thy1_idx = 1:7;
ai94_idx = 8:13;
camk_idx = 14:25;

for j = 1:13%length(data_list{1})
    data_dir = data_list{1}{j};
    [mouse_root_dir, exp_date, ~] = fileparts(data_dir);
    [~, mouse_id, ~] = fileparts(mouse_root_dir);
    [expt_root_dir, ~, ~] = fileparts(mouse_root_dir);
   
    grooming_file = [data_dir, filesep, 'grooming_events_roi_filtered.mat'];

    try
        disp('Loading grooming events')
        load(grooming_file)
    catch
        disp(['No grooming file for ', mouse_id, ' on ', exp_date, '...']);
        continue
    end
    groom_data{j} = all_events;

%         separate ipsilateral, contralateral, and bilateral events
%     ipsi{j} = zeros(size(FL_L));
%     contra{j} = zeros(size(FL_L));
%     bilat{j} = zeros(size(FL_L));
    ipsi(j) = 0;
    contra(j) = 0;
    bilat(j) = 0;
    for i = 1:size(event_type, 1)
        if contains(event_type(i,:), 'ipsi')
            ipsi(j) = ipsi(j)+1;
%             ipsi{j}(all_events(i,1):all_events(i,2)) = 1;
        elseif contains(event_type(i,:), 'cont')
            contra(j) = contra(j) + 1;
%             contra{j}(all_events(i,1):all_events(i,2)) = 1;
        elseif contains(event_type(i,:), 'bi')
            bilat(j) = bilat(j)+1;
%             bilat{j}(all_events(i,1):all_events(i,2)) = 1;
        else
            disp('wrong string in event type')
        end
    end


end
    
spon_idx = [1,2, 8, 9, 14, 15, 20, 21, 26, 27];
spon_idx = [1,2, 8, 9]; %test
evoked_idx = 1:length(groom_data);
evoked_idx(spon_idx) = [];



spon_groom = groom_data(spon_idx);
evoked_groom = groom_data(evoked_idx);
empty_idx = cellfun(@isempty, evoked_groom);
evoked_groom = evoked_groom(~empty_idx);


%%
spon_ipsi = ipsi(spon_idx);
evoked_ipsi = ipsi(evoked_idx);
evoked_ipsi = evoked_ipsi(~empty_idx);

%%
% evoked_ipsi = evoked_ipsi(~cellfun(@isempty, evoked_ipsi));

spon_contra = contra(spon_idx);
evoked_contra = contra(evoked_idx);
evoked_contra = evoked_contra(~empty_idx);
% evoked_contra = evoked_contra(~cellfun(@isempty, evoked_contra));

spon_bilat = bilat(spon_idx);
evoked_bilat = bilat(evoked_idx);
evoked_bilat = evoked_bilat(~empty_idx);
% evoked_bilat = evoked_bilat(~cellfun(@isempty, evoked_bilat));

% spon_ipsi = cellfun(@sum, spon_ipsi);
% spon_contra = cellfun(@sum, spon_contra);
% spon_bilat = cellfun(@sum, spon_bilat);
% 
% evoked_ipsi = cellfun(@sum, evoked_ipsi);
% evoked_contra = cellfun(@sum, evoked_contra);
% evoked_bilat = cellfun(@sum, evoked_bilat);

spon_all = [spon_ipsi; spon_contra; spon_bilat];
evoked_all = [evoked_ipsi; evoked_contra; evoked_bilat];

%%
max_trial_length = max([cellfun(@max, cellfun(@max, spon_groom, 'UniformOutput', false)), ...
    cellfun(@max, cellfun(@max, evoked_groom, 'UniformOutput', false))]);

event_raster = zeros(length(spon_groom) + length(evoked_groom), max_trial_length);
for i = 1:length(spon_groom)
    for j = 1:size(spon_groom{i},1)
        event_raster(i, spon_groom{i}(j,1):spon_groom{i}(j,2)) = 1;
    end
end
for i = 1:length(evoked_groom)
    for j = 1:size(evoked_groom{i}, 1)
        event_raster(i+length(spon_groom), evoked_groom{i}(j,1):evoked_groom{i}(j,2)) = 1;
    end
end

figure, 
t = xt(event_raster,fs,2);
subplot(3,3,1:2)
pcolor(t, size(event_raster,1):-1:1, 1-event_raster), 
colormap gray, shading flat
y = [size(event_raster,1), size(event_raster,1), ...
    size(event_raster,1)-length(spon_groom), size(event_raster,1)-length(spon_groom)];
patch([0 t(end) t(end) 0], y,'red', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
y = [size(event_raster,1)-length(spon_groom), size(event_raster,1)-length(spon_groom), 0, 0];
patch([0 t(end) t(end) 0], y,'blue', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
% line([0 0], ...
%     [size(event_raster, 1) size(event_raster,1)-length(spon_groom)], ...
%     'color', [1 0 0], ...
%     'LineWidth', 2);
% line([0 0], [0 length(evoked_groom)], 'color', [0 0 1], 'LineWidth', 2)
yticks([length(evoked_groom)/2, length(evoked_groom)+length(spon_groom)/2]);
yticklabels({'Evoked', 'Spontaneous'});
% xticks([])
% xlabel('Time (s)')
title('Grooming Events')
ax = gca;
ax.FontSize = 12;

num_spon_events = cellfun('size', spon_groom, 1);
num_evoked_events = cellfun('size', evoked_groom, 1);
[h,p] = ttest2(num_spon_events, num_evoked_events);


size_diff = length(num_evoked_events) - length(num_spon_events);
data2plot = [num_spon_events, nan(1, size_diff); num_evoked_events]';

subplot(3,3,[3,6]), 
boxplot(data2plot, 'colors', 'rb', 'Labels',{'Spontaneous','Evoked'});
ylabel('# Grooming events')
ax = gca;
ax.FontSize = 12;
title(['p=', num2str(p)])

% spon_duration = sum(event_raster(1:length(spon_groom),:), 2);
% evoked_duration = sum(event_raster(length(spon_groom)+1:end,:), 2);
% data2plot = [cat(1, spon_duration, nan(size_diff, 1)), evoked_duration]/fs;
% [h,p] = ttest2(spon_duration, evoked_duration);
% 
% 
% subplot(1,4,4), 
% boxplot(data2plot, 'colors', 'rb', 'Labels',{'Spontaneous','Evoked'});
% ylabel('Grooming duration (s)')
% ax = gca;
% ax.FontSize = 12;
% title(['p=', num2str(p)])


spon_probability = mean(event_raster(1:length(spon_groom), :));
spon_probability = movmean(spon_probability, round(5*fs));
evoked_probability = mean(event_raster(length(spon_groom)+1:end, :));
evoked_probability = movmean(evoked_probability, round(5*fs));
subplot(3,3,4:5)

% Get coefficients of a line fit through the data.
coefficients = polyfit(t, evoked_probability, 1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(t), max(t), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);
% Plot everything.
plot(xFit, yFit, 'k-', 'LineWidth', 1); % Plot fitted line.
hold on

plot(t, evoked_probability, 'b', 'LineWidth', 1.5)
plot(t, spon_probability, 'r', 'LineWidth', 1.5)
vline(32:60:t(end), 'b:')

axis tight
ylabel('Grooming probability')
xlabel('Time (s)')
ax = gca;
ax.FontSize = 12;


subplot(3,3,7)
names = ["Right", "Left", "Bilateral"];

p=pie(mean(spon_all, 2));
title('Spontaneous')
ax = gca;
ax.FontSize = 12;

subplot(3,3,8), pie(mean(evoked_all, 2)), title('Evoked')
ax = gca;
ax.FontSize = 12;
% subplot(3,3,8)
% lgd = legend(names);
