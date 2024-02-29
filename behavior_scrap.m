%% behavior analysis scrap
% 


clear, clc
fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
data_list = textscan(fileID, formatSpec);

fs = 90;

thy1_idx = 1:7;
ai94_idx = 8:13;
camk_idx = 14:25;

% create transition matrix
% states = ["Stop", "Elliptical", "Bilateral Asymmetric", "Bilateral", "Unilateral", "Lick"];
states = ["Stop", "Elliptical", "Asymmetric", "Bilateral", "Unilateral"];


spontaneous = [1,2,8,9,14,15,20,21];
evoked = [3:7,10:13,16:19,22:25];

%%

sessions = {'spontaneous', 'evoked'};

for i = 1:2
    tmat = zeros(numel(states));
    event_raster = [];

    for j = eval(sessions{i})%14:length(data_list{1})
        data_dir = data_list{1}{j};
        [mouse_root_dir, exp_date, ~] = fileparts(data_dir);
        [~, mouse_id, ~] = fileparts(mouse_root_dir);
        [expt_root_dir, ~, ~] = fileparts(mouse_root_dir);
       
        snippets_dir = [data_dir, filesep, 'snippets'];
        [snippets, labels] = parse_snippets(snippets_dir);
        num_events = cellfun(@(data) size(data, 1), snippets);
    
        event_duration{j} = cellfun(@(data) diff(data,1,2), snippets, 'UniformOutput', false);
    
        bmat = zeros(1, max(catcell(2, cellfun(@max, snippets, 'UniformOutput', false))));
        for ii = 1:length(snippets)
            for jj = 1:size(snippets{ii},1)
                bmat(snippets{ii}(jj,1):snippets{ii}(jj,2)) = ii;
            end
        end
    
        % consolidate unilateral events
        bmat2 = bmat;
        bmat2(bmat2==3)=2;
        bmat2(bmat2==4)=3;
        bmat2(bmat==5 | bmat==6) = 4;
        bmat2(bmat==7) = 0;% 5;
        bmat2 = bmat2 + 1;
        
        [episodes, idx] = aggregate(bmat, 3);
        

        for ii = 1:size(idx,1)
            tmp = bmat2(idx(ii,1):idx(ii,2));
            event_idx = [1 tmp] > 1;
            event_start = find(diff(event_idx) == 1);
            for jj = 1:length(event_start)
                if jj == 1
                    tmat(1, tmp(event_start(jj))) = tmat(1, tmp(event_start(jj)))+1;
                else
                    tmat(last_event, tmp(event_start(jj))) = tmat(last_event, tmp(event_start(jj))) + 1;
                end
                last_event = tmp(event_start(jj));
            end
            tmat(last_event, 1) = tmat(last_event, 1) + 1;
        end


    end



    figure
    subplot(1,2,1)
    imagesc(tmat)
    xticks(1:7), xticklabels(states)
    yticks(1:7), yticklabels(states)
    ylabel('From')
    xlabel('To')
    colormap(flipud(colormap('gray')))
    
    title('Transitions')
    
    colorbar;
    
    state_count = sum(tmat,2);
    pmat = tmat./state_count;
    
    subplot(1,2,2)
    imagesc(pmat)
    xticks(1:7), xticklabels(states)
    yticks(1:7), yticklabels(states)
    ylabel('From')
    xlabel('To')
    colormap(flipud(colormap('gray')))
    % clim([0 1])
    colorbar
    title('Conditional Probability')

    rel_freq(:,i) = state_count(2:end)./sum(state_count(2:end));

%     suptitle(sessions{i})
end



% % This is not working
% beta_num = tmat./sum(tmat(:));
% beta_denom = prod(state_count/sum(state_count));
% B = log(beta_num ./ beta_denom);
% 
% subplot(1,2,2)
% imagesc(B)
% xticks(1:7), xticklabels(states)
% yticks(1:7), yticklabels(states)
% ylabel('From')
% xlabel('To')
% colormap(flipud(colormap('gray')))
% % clim([0 1])
% colorbar
% title('Transition Bias')

% figure
% p = pie(state_count(2:end), states(2:end));
% fontsize(gcf, 12, "points")
% arrayfun(@(data) set(data, 'FaceAlpha', 0.5), p(1:2:end))



% sort rel_freq so unilateral
%%
figure, 
b = bar(rel_freq);
b(1).FaceColor = 'k';
b(2).FaceColor = [0.3010 0.7450 0.9330];
xticklabels(states(2:end))
ylabel('Relative Frequency')
legend({'Spontaneous', 'Evoked'}, 'Location', 'NorthEast', 'Box', 'Off')
set(gca, 'FontSize', 12)



%% find multiple events within snippets

event_duration_all = cell(1,7);
for i = 1:length(event_duration)
    if ~isempty(event_duration{i})
        for j = 1:length(event_duration{i})
            event_duration_all{j} = cat(1,event_duration_all{j}, event_duration{i}{j});
        end
    else
        disp('hello')
    end
end


figure, hold on
for i = 1:length(event_duration_all)
    subplot(4, 2, i)
    if any(event_duration_all{i}<0)
        event_duration_all{i} = event_duration_all{i}(event_duration_all{i}>0);
    end
    [f,xi] = ksdensity(event_duration_all{i}); 
    plot(xi./fs, f)
    findpeaks(f, fs, 'MinPeakProminence', 0.005)
%     histogram(event_duration_all{i}./fs, 50);
%     vline(median(event_duration_all{i})/fs, 'r')
%     vline(2*median(event_duration_all{i})/fs, 'r:')
%     vline(3*median(event_duration_all{i})/fs, 'r:')
%     vline(4*median(event_duration_all{i})/fs, 'r:')
    xlim([0 1])
    xlabel('Duration (s)')
    title(labels(i))
    ylabel('Counts')
end



%%
% check if events within grooming episodes follows a pattern
longest_event = max(diff(idx));
b = nan(size(idx,2), longest_event, 3);
for ii = 1:size(b,1)
    tmp = bmat(idx(1, ii): idx(2, ii));
    ellip = tmp == 1;
    unilat = tmp == 2 | tmp == 3;
    bilat = tmp == 4;

    b(ii, 1:length(ellip), 1) = ellip;
    b(ii, 1:length(unilat), 2) = unilat;
    b(ii, 1:length(bilat), 3) = bilat;
end





%% network visualization

G = digraph(pmat,states);
figure, %plot(G)
p=plot(G, 'LineWidth', G.Edges.Weight*10);
numlinks = state_count*1;
p.MarkerSize = numlinks;
p.NodeFontSize = 12;
p.ArrowSize = 20;
p.EdgeAlpha = 0.25;
%%

%         for ii = 1:size(idx,1)
%             tmp = bmat2(idx(ii, 1): idx(ii, 2));
%             for jj = 1:length(tmp)
%                 if jj == 1 && tmp(jj) ~= 0
%                     tmat(1, tmp(jj)) = tmat(1, tmp(jj))+1;
%                     last_event = tmp(jj);
%                 else
%                      if tmp(jj)~=tmp(jj-1) && tmp(jj) ~= 1
%                          tmat(last_event, tmp(jj)) = tmat(last_event, tmp(jj)) + 1;
%                          last_event = tmp(jj);
%                      end
%                 end
%             end
%             tmat(last_event, 1) = tmat(last_event, 1) + 1;
%         
%         end
    
    
    
    %     figure,   plot(xt(bmat, fs), bmat), yticks(1:7), yticklabels(labels)
    % 
    %     bmat2 = bmat;
    %     bmat2(bmat2==3)=2;
    %     bmat2(bmat2==4)=3;
    %     bmat2(bmat2>4) = 0;
    %     figure,   plot(xt(bmat, fs), bmat2), yticks(1:3), yticklabels(["elliptical", "unilateral", "bilateral"])
    
    





%%

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



