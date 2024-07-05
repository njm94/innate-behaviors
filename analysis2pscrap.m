clear, clc

% addpath(genpath('C:\Users\user\Documents\Nick\grooming'));

[event_file, event_path] = uigetfile('*.tsv','Select BORIS event labels.', 'Y:\nick\behavior\grooming\2p');
events = read_boris([event_path, filesep, event_file]);

[neuron_file, neuron_path] = uigetfile('*clean.mat','Select cleaned neuron data.', event_path);
load([neuron_path, filesep, neuron_file]);


%%
Nresample = resample(N, size(events,1), size(N, 2), 'Dimension', 2);


%%




% consolidate lick events
lick_idx = contains(events.Properties.VariableNames, 'Lick');
lick_events = events(:,lick_idx);
lick_events = any(table2array(lick_events),2);

% remove lick and point events from event matrix
idx = contains(events.Properties.VariableNames, 'Lick') | ...
    contains(events.Properties.VariableNames, 'Drop') | ...
    contains(events.Properties.VariableNames, 'Video');
stroke_events = removevars(events, idx);
%%

eventsarr = table2array(stroke_events)';

%%

figure, imagesc(imadjust([eventsarr; lick_events'; Nresample])),
colormap(flipud(colormap('gray')))



%%
fs = 90;
Ndenoise = smoothdata(Nresample, 2, 'lowess', [fs fs]);

clear r
tic
for j = 1:size(eventsarr,1)
    for i = 1:size(Nresample,1)
        r(j,i) = corr(eventsarr(j,:)', Ndenoise(i,:)');
    end
end
toc

%%

figure,
imagesc(r), 
colormap(redblue()), caxis([-0.3 0.3])
c=colorbar;
c.Label.String = 'Correlation';
yticklabels(stroke_events.Properties.VariableNames)
xlabel('Neuron Number')

%%
figure
for i = 1:size(r,1)
    for j = 1
        [rtmp, I] = max(r(i,:));
        if rtmp > 0.2
            idx = arr2idx(aggregate(eventsarr(i,:),3,fs));
            for k = 1:size(idx,1)
                subplot(1, size(r,1), i)
                plot(Ndenoise(I, idx(k,1)-450:idx(k,1)+450)), hold on
            end
        end
    end


end

%% plot location of highly correlated neurons onto the brain

close all
atlas_aligned_fig = fullfile(neuron_path, '..', 'dalsa', 'atlas_aligned.fig');

uiopen(atlas_aligned_fig,1)
hold on


for i = 1:size(N,1)
    x0 = double(cstat{i}.med(2));
    y0 = double(cstat{i}.med(1));

    [x,y] = get_neuron_location_in_dorsal_space( x0,y0, ...
        tforms{cstat{i}.use_tform}.roi_to_linear, ...
        tforms{cstat{i}.use_tform}.linear_to_wfield, ...
        tforms{cstat{i}.use_tform}.wfield_to_atlas);

    plot(x, y, 'k.', 'MarkerSize', 10)
end

prev_I = [];
tmp = r(6,:);
for i = 1:20
    [rr, I] = max(tmp);
    if rr<0.2, continue; end
%     plot(t, Nresample(I,:)-10*i, 'k') % "Elliptical" correlated neuron
    prev_I = cat(1, prev_I, I);
    tmp(I) = nan;
    x0 = double(cstat{I}.med(2));
    y0 = double(cstat{I}.med(1));

    [x,y] = get_neuron_location_in_dorsal_space( x0,y0, ...
        tforms{cstat{I}.use_tform}.roi_to_linear, ...
        tforms{cstat{I}.use_tform}.linear_to_wfield, ...
        tforms{cstat{I}.use_tform}.wfield_to_atlas);

    plot(x, y, 'r.', 'MarkerSize', 20)

end


tmp = r(1,:);
for i = 1:20
    [rr, I] = max(tmp);
    if rr<0.2, continue; end
%     plot(t, Nresample(I,:)-10*i, 'k') % "Elliptical" correlated neuron
    tmp(I) = nan;
    x0 = double(cstat{I}.med(2));
    y0 = double(cstat{I}.med(1));

    [x,y] = get_neuron_location_in_dorsal_space( x0,y0, ...
        tforms{cstat{I}.use_tform}.roi_to_linear, ...
        tforms{cstat{I}.use_tform}.linear_to_wfield, ...
        tforms{cstat{I}.use_tform}.wfield_to_atlas);
    
    if any(prev_I == I)
        plot(x,y, 'c.', 'MarkerSize', 20)
    else
        plot(x, y, 'b.', 'MarkerSize', 20)
    end
end


% axis tight
% yticklabels(fliplr(nlabel));
% title('Elliptical')






%% plot some example neurons

t = xt(Nresample, 90);
clc
% num_cells = 20l
nlabel = cell(1,10);
figure, subplot(1,2,1), hold on
tmp = r(1,:);
for i = 1:10
    [~, I] = max(tmp);
    plot(t, Nresample(I,:)-10*i, 'k') % "Elliptical" correlated neuron
    tmp(I) = nan;
    nlabel{1,i} = [nloc{I},'-N',num2str(I)];
    
end
axis tight
yticklabels(fliplr(nlabel));
title('Elliptical')

vline(t(find(table2array(events(:,3)))), 'r:')

patchplot(t(arr2idx(aggregate(eventsarr(1,:)', 3))), ylim, 'c', 0.5)
patchplot(t(arr2idx(aggregate(eventsarr(6,:)',3))), ylim, 'm', 0.5)
xlabel('Time (s)')


nlabel = cell(1,10);
subplot(1,2,2), hold on
tmp = r(6,:);
for i = 1:10
    [~, I] = max(tmp);
    plot(t, Nresample(I,:)-10*i, 'k') % "right" correlated neuron
    tmp(I) = nan;
    nlabel{1,i} = [nloc{I},'-N',num2str(I)];
end
axis([0 t(end) -105 0])

yticklabels(fliplr(nlabel));
title('Right')
vline(t(find(table2array(events(:,3)))), 'r:')

patchplot(t(arr2idx(aggregate(eventsarr(1,:)', 3))), ylim, 'c', 0.5)
patchplot(t(arr2idx(aggregate(eventsarr(6,:)',3))), ylim, 'm', 0.5)
xlabel('Time (s)')

%%
[~, I] = max(r(7,:))
plot(t, Nresample(I,:)-10, 'k') % "Right" correlated neuron

% [~, I] = max(r(8,:))
% plot(t, Nresample(I,:)-20, 'k') % "Left" correlated neuron


patchplot(t(arr2idx(aggregate(eventsarr(2,:), 3))), [-15 10], 'c', 0.5)
patchplot(t(arr2idx(aggregate(eventsarr(7,:),3))), [-15 10], 'm', 0.5)

xlabel('Time (s)')
yticks([])
% yticklabels(["Right", "Elliptical"])
ylabel('Example Neurons')
legend({'', '', 'Elliptical', 'Right'})
% patchplot(t(arr2idx(aggregate(eventsarr(8,:),3))), [-25 10], 'g', 0.5)


axis tight

%% neuron-behavior correlation matrix

stroke_types = table2array(events(:,3:end));
lick_consolidated = any(stroke_types(:,5:8), 2);
stroke_types(:,5:8) = [];
stroke_types = [stroke_types(:,1:4), lick_consolidated, stroke_types(:,5:end)];

%%
figure
thresh = 0.1:0.05:0.3;
for tt = 1:length(thresh)
    rmat = zeros(size(Ndenoise,1), size(stroke_types,2));
    for i = 1:size(Ndenoise,1)
        for j = 1:size(stroke_types,2)
            r_tmp = corr(Ndenoise(i,:)', stroke_types(:,j));
            if r_tmp > thresh(tt)
                rmat(i,j) = r_tmp;
            end
        end
    end
    
    
    %
    subplot(1,5,tt)
    test = sum(rmat>0, 2);
    histogram(test)
    ylabel('Number of neurons')
    if tt == 3
        xlabel('Number of behaviors neuron activity is highly correlated with')
    end
    title(['Thresh = ', num2str(thresh(tt))])
end



%%

figure, plot(Nresample(247,:))
