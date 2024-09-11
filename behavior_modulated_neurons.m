clear, clc, close all


% addpath(genpath('C:\Users\user\Documents\Nick\grooming'));
addpath('C:\Users\user\Documents\Nick\grooming\utils')



mp_list = {'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240731';
    'Y:\nick\behavior\grooming\2p\ECL3_thy1\20240802';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240729';
    'Y:\nick\behavior\grooming\2p\IDR3_tTA6s\20240731';
    'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240729';
    'Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240802'};


total_count = cell(1,length(mp_list));

for ii = 1:length(mp_list)
    [events, ~, event_table] = read_boris([mp_list{ii}, filesep, getAllFiles(mp_list{ii}, '.tsv')]);

    if ~isfile([mp_list{ii}, filesep, 'Nresample.mat'])
        load([mp_list{ii}, filesep, getAllFiles(mp_list{ii}, 'clean.mat')]);
        
        try Nresample = resample(N, size(events,1), size(N, 2), 'Dimension', 2);
        catch
            disp('Data too big. Splitting into halves and resampling each half separately')
            Nresample = resamplee(N', size(events,1), size(N,2))';
        end
    
        save([mp_list{ii}, 'Nresample.mat'], 'Nresample','nloc', 'cstat', 'tforms', '-v7.3')
    else
        disp('Loading resampled neuron data')
        load([mp_list{ii}, filesep, 'Nresample.mat'])
    end

    sz = 2;
    fs = 90;

    vel = readmatrix([mp_list{ii}, filesep, getAllFiles(mp_list{ii}, '_vel.csv')]);
    vid_end = find(events.("Video End"));
    
    vel = vel(1:vid_end,:);
    flrv = sum(vel(:,4:5).^2, 2).^0.5;
    fllv = sum(vel(:,7:8).^2, 2).^0.5;
    
    
    % consolidate lick events
    lick_idx = contains(events.Properties.VariableNames, 'Lick');
    lick_events = events(:,lick_idx);
    lick_events = any(table2array(lick_events),2);
    
    % remove lick and point events from event matrix
    idx = contains(events.Properties.VariableNames, 'Lick') | ...
        contains(events.Properties.VariableNames, 'Drop') | ...
        contains(events.Properties.VariableNames, 'Video');
    stroke_events = removevars(events, idx);
    behaviors = stroke_events.Properties.VariableNames;


    pos_mod = cell(1, length(behaviors));
    neg_mod = cell(1, length(behaviors));
    no_mod = cell(1,length(behaviors));
    
    pos_mod_loc = cell(1, length(behaviors));
    neg_mod_loc = cell(1, length(behaviors));
    no_mod_loc = cell(1, length(behaviors));
    
    
    for i = 1:length(behaviors)
        idx = strcmp(event_table.Behavior, behaviors{i}) & strcmp(event_table.BehaviorType, 'START');
        frame_start = event_table.ImageIndex(idx);
        clear tmp
    
    
    %     idx = arr2idx(aggregate(table2array(stroke_events(:,i)),1,fs));
    %     frame_start = idx(:,1);
    
        % This creates an array of size NxTxB,
        %   N (Number of neurons)
        %   T (Time centered on event start)
        %   B (Number of events)
        for j = 1:length(frame_start)
            try 
                tmp(:,:,j) = Nresample(:, frame_start(j)-(sz*fs):frame_start(j)+(sz*fs));
            catch
                disp('Inside catch')
            end
        end
        
        % Do a paired signrank test for each neuron's baseline vs response to 
        % determine if it is significantly modulated by the behavior
        no_mod_count = 1;
        neg_mod_count = 1;
        pos_mod_count = 1;
        pos_onset_count = 1;
        neg_onset_count = 1;
        for j = 1:size(tmp,1)
            baseline = squeeze(sum(tmp(j, 1:floor(size(tmp,2)/2), :),2));
            response = squeeze(sum(tmp(j, ceil(size(tmp,2)/2):end, :),2));
            onset_response = squeeze(sum(tmp(j, round(size(tmp,2)/2)-round(fs*sz/2):round(size(tmp,2)/2)+round(fs*sz/2), :),2));
            onset_baseline = squeeze(sum(tmp(j, [1:round(size(tmp,2)/2)-round(fs*sz/2), round(size(tmp,2)/2)+round(fs*sz/2):end], :),2));
            
    %         [p,h] = signrank(onset_response, onset_baseline);
    %         if h
    %             if mean(onset_baseline) > mean(onset_response)
    %                 neg_onset{i}(neg_onset_count,:,:) = tmp(j,:,:);
    %                 neg_onset_count = neg_onset_count + 1;
    %             else
    %                 pos_onset{i}(pos_onset_count,:,:) = tmp(j,:,:);
    %                 pos_onset_count = pos_onset_count + 1;
    %             end
    %         else
                [p, h] = signrank(baseline, response);
                if h
                    if mean(baseline) > mean(response)
                        neg_mod{i}(neg_mod_count,:,:) = tmp(j,:,:);
                        neg_mod_loc{i}(neg_mod_count) = nloc(j);
                        neg_mod_count = neg_mod_count + 1;
                    else
                        pos_mod{i}(pos_mod_count,:,:) = tmp(j,:,:);
                        pos_mod_loc{i}(pos_mod_count) = nloc(j);
                        pos_mod_count = pos_mod_count + 1;
                    end
                else
                    no_mod{i}(no_mod_count,:,:) = tmp(j,:,:);
                    if size(no_mod{i},3)>size(no_mod{i},2)
                        no_mod{i} = squeeze(no_mod{i})';
                    end
                    no_mod_loc{i}(no_mod_count) = nloc(j);
                    no_mod_count = no_mod_count + 1;
                    
                end
    %         end
            
        end
    
    end
    
    %%
    figure
    clear tmp
    
    caxlims = [-2 2];
    
    
    for i = 1:length(behaviors)
        tmp = mean(pos_mod{i},3);
        if ~isempty(tmp)    
            mr = mean(tmp(:,sz*fs:end), 2);
            
            [~, I] = sort(mr);   
            baseline = mean(tmp(:,1:sz*fs), 2);   
        
            subplot(6, length(behaviors), i)
            imagesc(tmp(I,:)-baseline(I))
    %         colorbar()
            clim(caxlims)
            colormap(bluewhitered())
            title(behaviors{i}), 
            xticks([])
    
            hold on
            vline(fs*sz, 'k-')
    
        
            subplot(6, length(behaviors), length(behaviors)+i)
            plot(xt(tmp, fs, 2)-sz, mean(tmp(I,:)-baseline(I)))
            hold on
            vline(0, 'k-')
        end
    %     sgsg
        
        tmp = mean(neg_mod{i},3);
    
        if ~isempty(tmp)    
            mr = mean(tmp(:,sz*fs:end), 2);
            [~, I] = sort(mr);
        
            baseline = mean(tmp(:,1:sz*fs), 2);   
        
            subplot(6, length(behaviors), 2*length(behaviors)+i)
            imagesc(tmp(I,:)-baseline(I))
            title(behaviors{i}), 
            xticks([])
    %         colorbar
            hold on
            vline(sz*fs, 'k-')
            clim(caxlims)
        
            subplot(6, length(behaviors), 3*length(behaviors)+i)
            plot(xt(tmp, fs, 2)-sz, mean(tmp(I,:)-baseline(I)))
            hold on
            vline(0, 'k-')
        end
    
        tmp = mean(no_mod{i},3);
    
        if ~isempty(tmp)    
            mr = mean(tmp(:,sz*fs:end), 2);
            [~, I] = sort(mr);
        
            baseline = mean(tmp(:,1:sz*fs), 2);   
        
            subplot(6, length(behaviors), 4*length(behaviors)+i)
            imagesc(tmp(I,:,:)-baseline(I))
            title(behaviors{i})
            clim(caxlims)
            hold on
            vline(sz*fs, 'k-')
            xticks([])
    
        
            subplot(6, length(behaviors), 5*length(behaviors)+i)
            plot(xt(tmp, fs, 2)-sz, mean(tmp(I,:)-baseline(I)))
            hold on
            vline(0, 'k-')
        end
    
    end
    caxis([-2 2])
    colormap(bluewhitered())

end

%%





%% plot single neuron responses to each behavior
close all
behaviors = stroke_events.Properties.VariableNames;


pos_mod = cell(1, length(behaviors));
neg_mod = cell(1, length(behaviors));
no_mod = cell(1,length(behaviors));

pos_mod_loc = cell(1, length(behaviors));
neg_mod_loc = cell(1, length(behaviors));
no_mod_loc = cell(1, length(behaviors));


for i = 1:length(behaviors)
    idx = strcmp(event_table.Behavior, behaviors{i}) & strcmp(event_table.BehaviorType, 'START');
    frame_start = event_table.ImageIndex(idx);
    clear tmp


%     idx = arr2idx(aggregate(table2array(stroke_events(:,i)),1,fs));
%     frame_start = idx(:,1);

    % This creates an array of size NxTxB,
    %   N (Number of neurons)
    %   T (Time centered on event start)
    %   B (Number of events)
    for j = 1:length(frame_start)
        try 
            tmp(:,:,j) = Nresample(:, frame_start(j)-(sz*fs):frame_start(j)+(sz*fs));
        catch
            disp('Inside catch')
        end
    end
    
    % Do a paired signrank test for each neuron's baseline vs response to 
    % determine if it is significantly modulated by the behavior
    no_mod_count = 1;
    neg_mod_count = 1;
    pos_mod_count = 1;
    pos_onset_count = 1;
    neg_onset_count = 1;
    for j = 1:size(tmp,1)
        baseline = squeeze(sum(tmp(j, 1:floor(size(tmp,2)/2), :),2));
        response = squeeze(sum(tmp(j, ceil(size(tmp,2)/2):end, :),2));
        onset_response = squeeze(sum(tmp(j, round(size(tmp,2)/2)-round(fs*sz/2):round(size(tmp,2)/2)+round(fs*sz/2), :),2));
        onset_baseline = squeeze(sum(tmp(j, [1:round(size(tmp,2)/2)-round(fs*sz/2), round(size(tmp,2)/2)+round(fs*sz/2):end], :),2));
        
%         [p,h] = signrank(onset_response, onset_baseline);
%         if h
%             if mean(onset_baseline) > mean(onset_response)
%                 neg_onset{i}(neg_onset_count,:,:) = tmp(j,:,:);
%                 neg_onset_count = neg_onset_count + 1;
%             else
%                 pos_onset{i}(pos_onset_count,:,:) = tmp(j,:,:);
%                 pos_onset_count = pos_onset_count + 1;
%             end
%         else
            [p, h] = signrank(baseline, response);
            if h
                if mean(baseline) > mean(response)
                    neg_mod{i}(neg_mod_count,:,:) = tmp(j,:,:);
                    neg_mod_loc{i}(neg_mod_count) = nloc(j);
                    neg_mod_count = neg_mod_count + 1;
                else
                    pos_mod{i}(pos_mod_count,:,:) = tmp(j,:,:);
                    pos_mod_loc{i}(pos_mod_count) = nloc(j);
                    pos_mod_count = pos_mod_count + 1;
                end
            else
                no_mod{i}(no_mod_count,:,:) = tmp(j,:,:);
                if size(no_mod{i},3)>size(no_mod{i},2)
                    no_mod{i} = squeeze(no_mod{i})';
                end
                no_mod_loc{i}(no_mod_count) = nloc(j);
                no_mod_count = no_mod_count + 1;
                
            end
%         end
        
    end

end

%%
figure
clear tmp

caxlims = [-2 2];


for i = 1:length(behaviors)
    tmp = mean(pos_mod{i},3);
    if ~isempty(tmp)    
        mr = mean(tmp(:,sz*fs:end), 2);
        
        [~, I] = sort(mr);   
        baseline = mean(tmp(:,1:sz*fs), 2);   
    
        subplot(6, length(behaviors), i)
        imagesc(tmp(I,:)-baseline(I))
%         colorbar()
        clim(caxlims)
        colormap(bluewhitered())
        title(behaviors{i}), 
        xticks([])

        hold on
        vline(fs*sz, 'k-')

    
        subplot(6, length(behaviors), length(behaviors)+i)
        plot(xt(tmp, fs, 2)-sz, mean(tmp(I,:)-baseline(I)))
        hold on
        vline(0, 'k-')
    end
%     sgsg
    
    tmp = mean(neg_mod{i},3);

    if ~isempty(tmp)    
        mr = mean(tmp(:,sz*fs:end), 2);
        [~, I] = sort(mr);
    
        baseline = mean(tmp(:,1:sz*fs), 2);   
    
        subplot(6, length(behaviors), 2*length(behaviors)+i)
        imagesc(tmp(I,:)-baseline(I))
        title(behaviors{i}), 
        xticks([])
%         colorbar
        hold on
        vline(sz*fs, 'k-')
        clim(caxlims)
    
        subplot(6, length(behaviors), 3*length(behaviors)+i)
        plot(xt(tmp, fs, 2)-sz, mean(tmp(I,:)-baseline(I)))
        hold on
        vline(0, 'k-')
    end

    tmp = mean(no_mod{i},3);

    if ~isempty(tmp)    
        mr = mean(tmp(:,sz*fs:end), 2);
        [~, I] = sort(mr);
    
        baseline = mean(tmp(:,1:sz*fs), 2);   
    
        subplot(6, length(behaviors), 4*length(behaviors)+i)
        imagesc(tmp(I,:,:)-baseline(I))
        title(behaviors{i})
        clim(caxlims)
        hold on
        vline(sz*fs, 'k-')
        xticks([])

    
        subplot(6, length(behaviors), 5*length(behaviors)+i)
        plot(xt(tmp, fs, 2)-sz, mean(tmp(I,:)-baseline(I)))
        hold on
        vline(0, 'k-')
    end

end
caxis([-2 2])
colormap(bluewhitered())




%% plot location of neurons on the brain

% close all
% clc
atlas_aligned_fig = fullfile(event_path, '..', 'dalsa', 'atlas_aligned.fig');

uiopen(atlas_aligned_fig,1)
hold on

for i = 1:size(Nresample,1)
    x0 = double(cstat{i}.med(2));
    y0 = double(cstat{i}.med(1));

    [x,y] = get_neuron_location_in_dorsal_space( x0,y0, ...
        tforms{cstat{i}.use_tform}.roi_to_linear, ...
        tforms{cstat{i}.use_tform}.linear_to_wfield, ...
        tforms{cstat{i}.use_tform}.wfield_to_atlas);

    plot(x, y, 'k.', 'MarkerSize', 5)
end

%% color atlas by modulation

load('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat');
load('C:\Users\user\Documents\Nick\grooming\utils\atlas.mat')

%%

figure
for i = 1:length(behaviors)
    tmp = zeros(size(dorsalMaps.dorsalMapScaled));
    mp = round(size(tmp,2)./2);
    if isempty(pos_mod_loc{i})
        continue
    end
    for j = 1:length(pos_mod_loc{i})    
        
        % If the value is negative, it is on the right side
        addmask = dorsalMaps.dorsalMapScaled == abs(areanames.(pos_mod_loc{i}{j}));
        if areanames.(pos_mod_loc{i}{j}) < 0
            addmask(:, 1:mp) = 0;
        else
            addmask(:,mp:end) = 0;
        end
        tmp = tmp + addmask;
    end

    [area, area_count] = unique(nloc);
    tmp2 = ones(size(dorsalMaps.dorsalMapScaled));
    for j = 1:length(area)
        % If the value is negative, it is on the right side
        addmask = (area_count(j)-1).*(dorsalMaps.dorsalMapScaled == abs(areanames.(area{j})));
        if areanames.(area{j}) < 0
            addmask(:, 1:mp) = 0;
        else
            addmask(:,mp:end) = 0;
        end
        tmp2 = tmp2 + addmask;

    end
    subplot(1,length(behaviors), i)
    imagesc(tmp./tmp2)
    title(behaviors{i})
    colorbar,
    xticks([]), yticks([])
    caxis([0 20])
    colormap(bluewhitered())
    hold on
    for p = 1:length(dorsalMaps.edgeOutline)
        plot(dorsalMaps.edgeOutline{p}(:,2), dorsalMaps.edgeOutline{p}(:,1), 'k')
    end

end


%% negatively modulated neurons
figure
for i = 1:length(behaviors)
    tmp = zeros(size(dorsalMaps.dorsalMapScaled));
    mp = round(size(tmp,2)./2);
    if isempty(neg_mod_loc{i})
        continue
    end
    for j = 1:length(neg_mod_loc{i})
        % If the value is negative, it is on the right side
        addmask = dorsalMaps.dorsalMapScaled == abs(areanames.(neg_mod_loc{i}{j}));
        if areanames.(neg_mod_loc{i}{j}) < 0
            addmask(:, 1:mp) = 0;
        else
            addmask(:,mp:end) = 0;
        end
        tmp = tmp - addmask;
    end

    [area, area_count] = unique(nloc);
    tmp2 = ones(size(dorsalMaps.dorsalMapScaled));
    for j = 1:length(area)
        % If the value is negative, it is on the right side
        addmask = (area_count(j)-1).*(dorsalMaps.dorsalMapScaled == abs(areanames.(area{j})));
        if areanames.(area{j}) < 0
            addmask(:, 1:mp) = 0;
        else
            addmask(:,mp:end) = 0;
        end
        tmp2 = tmp2 + addmask;

    end

    subplot(1,length(behaviors), i)
    imagesc(tmp./tmp2)
    title(behaviors{i})
    colorbar,
    xticks([]), yticks([])
    caxis([-20 0])
    colormap(bluewhitered())
    hold on
    for p = 1:length(dorsalMaps.edgeOutline)
        plot(dorsalMaps.edgeOutline{p}(:,2), dorsalMaps.edgeOutline{p}(:,1), 'k')
    end

end

