% Regress movement out of 2p data

clear, clc

% addpath(genpath('C:\Users\user\Documents\Nick\grooming'));
addpath('C:\Users\user\Documents\Nick\grooming\utils')
[event_file, event_path] = uigetfile('*.tsv','Select BORIS event labels.', 'Y:\nick\behavior\grooming\2p');
[events, b_idx, ~] = read_boris([event_path, filesep, event_file]);

if ~isfile([event_path, 'Nresample.mat'])
    [neuron_file, neuron_path] = uigetfile('*clean.mat','Select cleaned neuron data.', event_path);
    load([neuron_path, filesep, neuron_file]);
    
    disp("Resampling neural data to match behavior")
    try Nresample = resample(N, size(events,1), size(N, 2), 'Dimension', 2);
    catch
        disp('Data too big. Splitting into halves and resampling each half separately')
        Nresample = resamplee(N', size(events,1), size(N,2))';
    end

    save([neuron_path, 'Nresample.mat'], 'Nresample', 'cstat', 'nloc', 'tforms')
else
    disp('Loading resampled neuron data')
    load([event_path, 'Nresample.mat'])
end


%% load DLC tracks
vel = readmatrix([event_path, filesep, getAllFiles(event_path, '_vel.csv')]);
vid_end = find(events.("Video End"));

vel = vel(1:vid_end,:);
flrv = sum(vel(:,4:5).^2, 2).^0.5;
fllv = sum(vel(:,7:8).^2, 2).^0.5;


%% load audio and valve open information

test = readmatrix([event_path, filesep, getAllFiles(event_path, '_trim.txt')]);
audio = test(1:vid_end,5);
valve = test(1:vid_end,6);

%%
flrthresh = flrv>mean(flrv) + std(flrv);
fllthresh = fllv>mean(fllv) + std(fllv);


clear Bmean
count = 1;
labs = {};

for i = 1:size(events,2)
    if contains(events.Properties.VariableNames{i}, 'Drop') || contains(events.Properties.VariableNames{i}, 'Video') || contains(events.Properties.VariableNames{i}, 'Lick')
        continue
    end
    event_table(:, count) = table2array(events(:,i));
    Bmean(:,count) = mean(Nresample(:,logical(event_table(:,count))),2);
    labs{count} = events.Properties.VariableNames{i};
    count = count + 1;
    
end

flrthresh(any(event_table, 2)) = 0;
fllthresh(any(event_table, 2)) = 0;


Bmean = cat(2, Bmean, mean(Nresample(:, logical(flrthresh)),2));
Bmean = cat(2, Bmean, mean(Nresample(:, logical(fllthresh)),2));
% Bmean = cat(2, Bmean, mean(Nresample(:, logical(audio)),2));
labs = [labs, 'FLR', 'FLL'];
% figure, imagesc(Bmean(I,:)), xticklabels(labs)
% colorbar
% caxis([-1 3])
% colormap(bluewhitered())


%% hierarchical clustering looks good



% distMatrix = pdist(Bmean', 'euclidean');
% Z = linkage(distMatrix, 'average');
Z = linkage(Bmean', 'average', 'correlation');
% Z = linkage(Bmean', 'average', 'correlation');

figure, subplot(3,2,1)
cutoff = 0.3;
[H, T, outperm] = dendrogram(Z,'Labels', labs);  % Plot the dendrogram
% [H, T, outperm] = dendrogram(Z,'Labels', labs, 'ColorThreshold', cutoff);  % Plot the dendrogram
ylabel('Distance (1-correlation)')
xticks([])


% sort neurons by anatomical location
label_column = ones(size(Bmean,1), 2);
anat = sort(unique(nloc));
I = [];
for i = 1:length(anat)
    Itmp = find(strcmp(nloc,anat{i}));
%     [~, Itmp_sorted] = sort(mean(Bmean(Itmp,end-1:end),2)); % forelimbs
    [~, Itmp_sorted] = sort(mean(Bmean(Itmp,[4,6]),2));
%     I = [I; find(strcmp(nloc, anat{i}))];
    I = [I; Itmp(Itmp_sorted)];

end


anat_parent = unique(cellfun(@(x) x(1:3), anat, 'UniformOutput', false));
for i = 1:length(anat_parent)
    label_column(contains(nloc(I), anat_parent{i}), 1) = i;  
end

for i = 1:length(anat)
    label_column(contains(nloc(I), anat{i}), 2) = i;    
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
subplot(3,2,[4, 6]), imagesc(label_column), colormap default


%%
% Determine the number of clusters (for example, 3 clusters)
numClusters = 2;

% Assign each observation to a cluster
clusterIdx = cluster(Z, 'maxclust', numClusters);

% Calculate silhouette values
silhouetteValues = silhouette(Bmean', clusterIdx);

% Plot the silhouette values
figure;
silhouette(Bmean', clusterIdx);
title('Silhouette Plot for Hierarchical Clustering');
xlabel('Silhouette Value');
ylabel('Cluster');


% Calculate the mean silhouette value for each cluster
meanSilhouette = zeros(numClusters, 1);
for i = 1:numClusters
    meanSilhouette(i) = mean(silhouetteValues(clusterIdx == i));
end

% Display the mean silhouette values for each cluster
disp('Mean Silhouette Values for each cluster:');
disp(meanSilhouette);

% Perform ANOVA to compare silhouette values across clusters
[p, tbl, stats] = anova1(silhouetteValues, clusterIdx);

% If p-value is significant, it indicates a difference in clustering quality
if p < 0.05
    disp('There is a significant difference in silhouette values between clusters.');
else
    disp('No significant difference in silhouette values between clusters.');
end

% Post-hoc test to determine which clusters differ
[c, m] = multcompare(stats);
%%

numBootstrap = 1000; % Number of bootstrap samples
bootstrapClusters = cell(numBootstrap, 1);

for i = 1:numBootstrap
    % Generate a bootstrap sample
    bootstrapSample = datasample(Bmean, size(Bmean, 1));
    
    % Perform hierarchical clustering on the bootstrap sample
    Z_bootstrap = linkage(bootstrapSample', 'ward');
    
    % Store the bootstrap clustering result
    bootstrapClusters{i} = Z_bootstrap;
end

%% 
% Initialize AU and BP values
AU = zeros(size(Z, 1), 1);
BP = zeros(size(Z, 1), 1);

% Iterate over the original clusters to calculate AU and BP
for j = 1:size(Z, 1)
    originalCluster = cluster(Z, 'maxclust', j);
    
    % Calculate BP as the proportion of bootstrap clusters that match the original cluster
    matchingClusters = 0;
    for i = 1:numBootstrap
        bootstrapCluster = cluster(bootstrapClusters{i}, 'maxclust', j);
        if isequal(bootstrapCluster, originalCluster)
            matchingClusters = matchingClusters + 1;
        end
    end
    
    BP(j) = matchingClusters / numBootstrap;
    
    % Approximate AU by adjusting BP with a method such as bootstrapping
    AU(j) = 2 * BP(j); % Simplified approximation (in practice, more sophisticated adjustments are made)
end



figure;
dendrogram(Z);
hold on;

% Annotate the dendrogram with AU and BP values
for j = 1:length(AU)
    text(j, Z(j,3), sprintf('AU: %.2f, BP: %.2f', AU(j), BP(j)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
hold off;
