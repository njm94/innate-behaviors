% Analyze changes in functional connectivity associated with long duration
% continuous grooming behaviors
%
% Expected outcome - Increase in FC at onset, followed by reduced FC at the
% end
%
% This is on 1p data only


%% Look at power spectra
fs = 90;
count = 1;
for i =1:length(roi_signal)
    if ~isempty(roi_signal{i})
        for j = 1:length(roi_signal{i})
            if ~isempty(roi_signal{i}{j})
                p_before(count,:) = 10*log10(mean(pspectrum(roi_signal{i}{j}(1:blen*fs,:), fs),2));
                p_early(count,:) = 10*log10(mean(pspectrum(roi_signal{i}{j}(1+blen*fs:2*blen*fs,:), fs),2));
                p_late(count,:) = 10*log10(mean(pspectrum(roi_signal{i}{j}(end-(2*blen*fs)+1:end-(blen*fs),:), fs),2));
                p_after(count,:) = 10*log10(mean(pspectrum(roi_signal{i}{j}(end-(blen*fs)+1:end,:), fs),2));
                count = count +1;
            end
        end
    end
end
[~, f] = pspectrum(roi_signal{i}{j}(1:blen*fs,:), fs);

%%
figure, hold on
%plot(f, mean(p_before))
plot(f, mean(p_early-p_before))
plot(f, mean(p_early-p_before)+std(p_early-p_before), ':')
plot(f, mean(p_late-p_before))
plot(f, mean(p_after-p_before))
%%

midx = {'thy1_idx', 'ai94_idx', 'hyl3_idx', 'ibl2_idx'};



thy1_idx = 1:12;
ai94_idx = 13:16;
hyl3_idx = 17:40;
ibl2_idx = 41:55;
% m = 4


test_pre = pre_corrmat;
% test_pre = test_pre(:,eval(midx{m}));
test_pre = catcell(3, test_pre(~cellfun(@isempty, test_pre)));



upper_t = triu(true(size(test_pre,1)), 1); % Mask for upper triangle
test_pre = reshape(test_pre(repmat(upper_t, 1, 1, size(test_pre,3))), [], size(test_pre,3));
test_pre = tanh(mean(atanh(test_pre)))';
% test_pre = mean(zscore(test_pre,[],2))';
% test_pre = mean(test_pre)';

test_early = early_corrmat;
% test_early = test_early(:,eval(midx{m}));
test_early = catcell(3, test_early(~cellfun(@isempty, test_early)));
test_early = reshape(test_early(repmat(upper_t, 1, 1, size(test_early,3))), [], size(test_early,3));
test_early = tanh(mean(atanh(test_early)))';
% test_early = mean(zscore(test_early,[],2))';
% test_early = mean(test_early)';

test_late = late_corrmat;
% test_late = test_late(:,eval(midx{m}));
test_late = catcell(3, test_late(~cellfun(@isempty, test_late)));
test_late = reshape(test_late(repmat(upper_t, 1, 1, size(test_late,3))), [], size(test_late,3));
test_late = tanh(mean(atanh(test_late)))';
% test_late = mean(zscore(test_late, [], 2))';
% test_late = mean(test_late)';


test_post = post_corrmat;
% test_post = test_post(:,eval(midx{m}));
test_post = catcell(3, test_post(~cellfun(@isempty, test_post)));
test_post = reshape(test_post(repmat(upper_t, 1, 1, size(test_post,3))), [], size(test_post,3));
test_post = tanh(mean(atanh(test_post)))';
% test_post = mean(zscore(test_post,[],2))';
% test_post = mean(test_post)';

% remove zeros for now but check where they are coming from in the future

if any(test_pre==0) 
    zidx=find(test_pre==0); 
    test_pre(zidx) = [];
    test_early(zidx) = [];
    test_late(zidx) = [];
    test_post(zidx) = [];
end

plot_data = [test_pre, test_early, test_late, test_post];

figure, boxplot(plot_data, 'Colors', 'k', 'Symbol', '')
hold on
swarmchart(repmat([1 2 3 4], size(test_early,1), 1), plot_data, 'k', 'XJitterWidth', 0.25)
cols = [0 0.4470 0.7410; 0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560];
for i = 1:length(midx)
    plot(mean(plot_data(eval(midx{i}),:))', 'Color', cols(i,:), 'LineWidth', 2)
end

ylabel('Average correlation')
xticklabels({'Before 5s', 'First 5s', 'Last 5s', 'After 5s'})
ax = gca;
ax.FontSize = 14;




%% better stats

% Convert data into long format
num_trials = size(plot_data, 1);
timepoints = repmat((1:4), num_trials, 1);  % Timepoint variable
mouseID = repmat([ones(1,12), 1+ones(1,4), 2+ones(1,24),3+ones(1,15)]', 1, size(timepoints,2));  % Mouse identity variable
response = plot_data(:);  % Flatten response matrix

% Rank-transform response for non-parametric analysis
ranked_response = tiedrank(response);

% Convert to table
T = table(mouseID(:), timepoints(:), ranked_response, 'VariableNames', {'Mouse', 'Timepoint', 'Response'});

% Fit LMM (random intercept for Mouse)
lme = fitlme(T, 'Response ~ Timepoint + (1|Mouse)', 'FitMethod', 'REML');

% Display results
disp(anova(lme));

%%
% Get estimated fixed effects and covariance matrix
beta = fixedEffects(lme); % Fixed effects (Intercept + Timepoint)
C = lme.CoefficientCovariance; % Covariance matrix

% Define contrast matrix for pairwise comparisons
% Each row represents a contrast: e.g., [0 1 -1 0] tests Timepoint 1 vs 2
contrasts = [0  1 -1  0;  % Timepoint 1 vs Timepoint 2
             0  1  0 -1;  % Timepoint 1 vs Timepoint 3
             0  0  1 -1;  % Timepoint 2 vs Timepoint 3
             0  1 -0 -1]; % Timepoint 1 vs Timepoint 4

% Compute p-values for each contrast
for i = 1:size(contrasts,1)
    pval = coefTest(lme, contrasts(i,:));
    fprintf('Comparison %d: p = %.4f\n', i, pval);
end
%%

% Run Friedman test
[p_friedman, tbl, stats] = friedman(plot_data, 1); % 1 indicates repeated measures
disp(['Friedman test p-value: ', num2str(p_friedman)]);

% If the Friedman test is significant, perform post-hoc pairwise comparisons
if p_friedman < 0.05
    fprintf('Post-hoc pairwise comparisons:\n');
    
    % Define time points
    time_labels = {'Pre', 'Early', 'Late', 'Post'};
    
    % Get all possible pairwise combinations
    pairs = nchoosek(1:4, 2);
    p_values = zeros(size(pairs, 1), 1);
    
    % Perform Wilcoxon signed-rank tests for each pair
    for i = 1:size(pairs, 1)
        [p_values(i), ~, stats] = signrank(plot_data(:, pairs(i, 1)), plot_data(:, pairs(i, 2)));
        fprintf('%s vs. %s: p = %.4f\n', time_labels{pairs(i, 1)}, time_labels{pairs(i, 2)}, p_values(i));
    end

    % Multiple comparisons correction (Bonferroni)
    corrected_p_values = min(p_values * length(p_values), 1); % Bonferroni correction
    fprintf('\nBonferroni corrected p-values:\n');
    
    for i = 1:size(pairs, 1)
        fprintf('%s vs. %s: p_corrected = %.4f\n', time_labels{pairs(i, 1)}, time_labels{pairs(i, 2)}, corrected_p_values(i));
    end
end

timepoints = {'test_pre', 'test_early', 'test_late', 'test_post'};
figure,  hold on
for i = 1:length(timepoints)
    [a,b] = ksdensity(eval(timepoints{i}));
    plot(a,b, 'LineWidth', 2);
end
legend({'Before 5s', 'First 5s', 'Last 5s', 'After 5s'}, 'Location', 'Best')



%%
% hyl3_idx = 14:19;
% ibl2_idx = 20:25;
midx = {'thy1_idx', 'ai94_idx', 'hyl3_idx', 'ibl2_idx'};
thy1_idx = 1:7;
ai94_idx = 8:13;
hyl3_idx = 14:19;
ibl2_idx = 20:25;

m = 3;
caxlims = [0.5 1];
test_pre = pre_corrmat;
test_pre = pre_corrmat(:,eval(midx{m}));
test_pre = catcell(3, test_pre(~cellfun(@isempty, test_pre)));
disp(size(test_pre,3))

test_early = early_corrmat;
test_early = early_corrmat(:,eval(midx{m}));
test_early = catcell(3, test_early(~cellfun(@isempty, test_early)));

test_late = late_corrmat;
test_late = late_corrmat(:,eval(midx{m}));
test_late = catcell(3, test_late(~cellfun(@isempty, test_late)));

test_post = post_corrmat;
test_post = post_corrmat(:,eval(midx{m}));
test_post = catcell(3, test_post(~cellfun(@isempty, test_post)));


figure
subplot(1,4,1), imagesc(tanh(nanmean(atanh(test_pre),3))), colorbar, caxis(caxlims)
xticks(1:size(test_pre,1))
yticks(1:size(test_pre,1))
xticklabels(labels)
yticklabels(labels)
subplot(1,4,2), imagesc(tanh(nanmean(atanh(test_early), 3))), colorbar, caxis(caxlims)
xticks(1:size(test_pre,1))
yticks(1:size(test_pre,1))
xticklabels(labels)
yticklabels(labels)
subplot(1,4,3), imagesc(tanh(nanmean(atanh(test_late), 3))), colorbar, caxis(caxlims)
xticks(1:size(test_pre,1))
yticks(1:size(test_pre,1))
xticklabels(labels)
yticklabels(labels)

subplot(1,4,4), imagesc(tanh(nanmean(atanh(test_post), 3))), colorbar, caxis(caxlims)
xticks(1:size(test_pre,1))
yticks(1:size(test_pre,1))
xticklabels(labels)
yticklabels(labels)

colormap(flipud(colormap(gray)))
%%

tb = tanh(mean(atanh(test_pre),3));
tb = tb(upper_t);

te = tanh(mean(atanh(test_early),3));
te = te(upper_t);

tl = tanh(mean(atanh(test_late),3));
tl = tl(upper_t);

ta = tanh(mean(atanh(test_post),3));
ta = ta(upper_t);


figure, plot(sort(tb), (1:numel(tb))/numel(tb), 'LineWidth', 2);
hold on
plot(sort(te), (1:numel(te))/numel(te), 'LineWidth', 2)
plot(sort(tl), (1:numel(tl))/numel(tl), 'LineWidth', 2)
plot(sort(ta), (1:numel(ta))/numel(ta), 'LineWidth', 2)
legend({'Before 5', 'First 5', 'Last 5', 'After 5'}, 'FontSize', 16, 'Location', 'Best')



%%
figure
d_early_pre = tanh(nanmean(atanh(test_early),3))-tanh(nanmean(atanh(test_pre), 3));
subplot(1,3,1), imagesc(d_early_pre), 
caxis([-0.2 0.2])
xticks(1:size(test_pre,1))
yticks(1:size(test_pre,1))
xticklabels(labels)
yticklabels(labels)
title('First 5s - Before 5s')
c=colorbar;
c.Label.String = '\Delta Correlation';

d_late_pre = tanh(nanmean(atanh(test_late),3))-tanh(nanmean(atanh(test_pre), 3));
subplot(1,3,2), imagesc(d_late_pre),  caxis([-0.2 0.2])
xticks(1:size(test_pre,1))
yticks(1:size(test_pre,1))
xticklabels(labels)
yticklabels(labels)
title('Last 5s - Before 5s')
c=colorbar;
c.Label.String = '\Delta Correlation';

d_post_pre = tanh(nanmean(atanh(test_post), 3))-tanh(nanmean(atanh(test_pre), 3));
subplot(1,3,3), imagesc(d_post_pre), colorbar, caxis([-0.2 0.2])
xticks(1:size(test_pre,1))
yticks(1:size(test_pre,1))
xticklabels(labels)
yticklabels(labels)
title('After 5s - Before 5s')
c=colorbar;
c.Label.String = '\Delta Correlation';

colormap(bluewhitered())
% for i = 1:length(data_list{1})
%     subplot(1,3,1)
%     imagesc()
% end



d_time_points = {'d_early_pre', 'd_late_pre', 'd_post_pre'};
time_points = {'test_early', 'test_late', 'test_post'};
% figure,

data_1 = atanh(test_pre);
% data_1 = test_pre;



for i = 1:length(time_points)
    d_tmp = eval(d_time_points{i});
    % data_2 = eval(time_points{i});
    data_2 = atanh(eval(time_points{i}));

    clear p_values
    % for ii = 1:size(data_1,1)
    %     for jj = 1:size(data_2,1)
    %         [p_values(ii,jj), observeddifference, effectsize] = permutationTest(squeeze(data_1(ii,jj,:)), squeeze(data_2(ii,jj,:)), 1000); 
    %     end
    % end

    for ii = 1:size(data_1,1)
        for jj = 1:size(data_2,1)
            [~, p_values(ii,jj)] = ttest(squeeze(data_1(ii,jj,:)), squeeze(data_2(ii,jj,:))); 
        end
    end


    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_values(upper_t),0.05,'dep','yes');
    h_tmp = zeros(size(p_values));
    h_tmp(upper_t) = h;

    if nansum(h(:)) ==0, h_tmp = zeros(size(d_tmp)); end
    G = graph(d_tmp.*h_tmp, 'upper');
    LWidths = abs(G.Edges.Weight)*20;
    LColors = zeros(numedges(G), 3);

    for ii = 1:numedges(G)
        if G.Edges.Weight(ii) <0
            LColors(ii,:) = [0 0 1];
        else
            LColors(ii,:) = [1 0 0];
        end
    end
    figure, hold on
    
    for p = 1:length(dorsalMaps.edgeOutline)
        plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k', 'LineWidth', 1);
        xticks([])
        yticks([])
    end
    set(gca, 'YDir', 'reverse');
    hold on, axis equal off
    h=plot(G, 'XData', seeds(:,1), 'YData', seeds(:,2), 'EdgeColor', LColors, 'LineWidth', LWidths, 'NodeLabel', {}, 'NodeColor', 'k');
end
%%


timepoints = {'test_early', 'test_late', 'test_post'};


for i = 1:length(timepoints)
    data_1 = atanh(test_pre);
    nanmask = eye(size(data_1,1), size(data_1,2));
    nanmask = repmat(nanmask, 1, 1, size(data_1, 3));
    nanmask(nanmask==1) = nan;

    data_1 = data_1 + nanmask;
    data_1 = tanh(squeeze(mean(data_1, 'omitnan')));

    data_2 = atanh(eval(timepoints{i}));
    data_2 = data_2 + nanmask;
    data_2 = tanh(squeeze(mean(data_2, 'omitnan')));

    data_1(:,31) = []; data_2(:,31) = []; % empty column

    tmp = data_2 - data_1;

    figure, bar(mean(tmp(:,42:end),2)), hold on
    % errorbar(mean(tmp(:,1:12),2),std(tmp(:,1:12),[],2), 'k.')
    grid on
    xticks(1:size(data_2,1));
    xticklabels(labels)
    axis([-0.2000 29.2000, -0.4 0.2])

    % figure, imagesc(data_2 - data_1)
    % ylabel('region')
    % xlabel('trial')
    % caxis([-0.2 0.2])
    % colormap(bluewhitered()),
    % colorbar
    
    % continue
    % sgfgsf
    % data_2 = eval(time_points{i});
    

    clear p_values
    for ii = 1:size(data_1,1)
        for jj = 1:size(data_2,1)
            % [p_values(ii,jj), observeddifference, effectsize] = permutationTest(squeeze(data_1(ii,jj,:)), squeeze(data_2(ii,jj,:)), 1000); 
            [p_values(ii,jj), observeddifference, effectsize] = permutationTest(squeeze(data_1(ii,:)), squeeze(data_2(jj,:)), 1000); 
        end
    end

    % for ii = 1:size(data_1,1)
    %     for jj = 1:size(data_2,1)
    %         % [~, p_values(ii,jj)] = ttest(squeeze(data_1(ii,jj,:)), squeeze(data_2(ii,jj,:))); 
    %         [~, p_values(ii,jj)] = ttest(squeeze(data_1(ii,:)), squeeze(data_2(jj,:))); 
    %     end
    % end


    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_values(upper_t),0.05,'dep','yes');
    h_tmp = zeros(size(p_values));
    h_tmp(upper_t) = h;

    if nansum(h(:)) ==0, h_tmp = zeros(size(d_tmp)); end
    G = graph(d_tmp.*h_tmp, 'upper');
    LWidths = abs(G.Edges.Weight)*20;
    LColors = zeros(numedges(G), 3);

    for ii = 1:numedges(G)
        if G.Edges.Weight(ii) <0
            LColors(ii,:) = [0 0 1];
        else
            LColors(ii,:) = [1 0 0];
        end
    end
    figure, hold on
    
    for p = 1:length(dorsalMaps.edgeOutline)
        plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k', 'LineWidth', 1);
        xticks([])
        yticks([])
    end
    set(gca, 'YDir', 'reverse');
    hold on, axis equal off
    h=plot(G, 'XData', seeds(:,1), 'YData', seeds(:,2), 'EdgeColor', LColors, 'LineWidth', LWidths, 'NodeLabel', {}, 'NodeColor', 'k');
end

