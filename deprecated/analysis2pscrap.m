% attempt to do parcellation on 1p signals            
tmp = dFF_crop(:,:,prewin);
tmp = reshape(tmp, size(tmp,1).*size(tmp,2), size(tmp,3));
k=5;
[kidx, ~] = kmeans(tmp, k, 'Distance', 'correlation');
kmap_pre{j}(:,:,i) = reshape(kidx, 128, 128);

tmp = dFF_crop(:,:,earlywin);
tmp = reshape(tmp, size(tmp,1).*size(tmp,2), size(tmp,3));
[kidx, ~] = kmeans(tmp, k, 'Distance', 'correlation');
kmap_early{j}(:,:,i) = reshape(kidx, 128, 128);

tmp = dFF_crop(:,:,latewin);
tmp = reshape(tmp, size(tmp,1).*size(tmp,2), size(tmp,3));
[kidx, ~] = kmeans(tmp, k, 'Distance', 'correlation');
kmap_late{j}(:,:,i) = reshape(kidx, 128, 128);

tmp = dFF_crop(:,:,postwin);
tmp = reshape(tmp, size(tmp,1).*size(tmp,2), size(tmp,3));
[kidx, ~] = kmeans(tmp, k, 'Distance', 'correlation');
kmap_post{j}(:,:,i) = reshape(kidx, 128, 128);
disp('Done with kmaps')
% continue



clear, clc

% addpath(genpath('C:\Users\user\Documents\Nick\grooming'));
addpath('C:\Users\user\Documents\Nick\grooming\utils')
[event_file, event_path] = uigetfile('*.tsv','Select BORIS event labels.', 'Y:\nick\behavior\grooming\2p');
[events, ~, event_table] = read_boris([event_path, filesep, event_file]);

if ~isfile([event_path, 'Nresample.mat'])
    [neuron_file, neuron_path] = uigetfile('*clean.mat','Select cleaned neuron data.', event_path);
    load([neuron_path, filesep, neuron_file]);
    
    
    try Nresample = resample(N, size(events,1), size(N, 2), 'Dimension', 2);
    catch
        disp('Data too big. Splitting into halves and resampling each half separately')
        Nresample = resamplee(N', size(events,1), size(N,2))';
    end

    save([neuron_path, 'Nresample.mat'], 'Nresample')
else
    disp('Loading resampled neuron data')
    load([event_path, 'Nresample.mat'])
end


vel = readmatrix([event_path, filesep, getAllFiles(event_path, '_vel.csv')]);
vid_end = find(events.("Video End"));

vel = vel(1:vid_end,:);
flrv = sum(vel(:,4:5).^2, 2).^0.5;
fllv = sum(vel(:,7:8).^2, 2).^0.5;
%%
% 
% 
% tmpf = U * Vm;
% fullMat = corr(Nresample', tmpf');
% 
% tmpr = U * VmR;
% redMat = corr(Nresample', tmpr');



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

figure, imagesc(imadjust([Nresample])),
colormap(flipud(colormap('gray')))



%%
fs = 90;
% Ndenoise = smoothdata(Nresample, 2, 'lowess', [fs fs]);

clear r
tic
for j = 1:size(eventsarr,1)
    for i = 1:size(Nresample,1)
        r(j,i) = corr(eventsarr(j,:)', Nresample(i,:)');
    end
end
toc

%% plot single neuron responses to each behavior
% close all
behaviors = stroke_events.Properties.VariableNames;
sz = 5;
figure
pos_mod = cell(1, length(behaviors));
neg_mod = cell(1, length(behaviors));
no_mod = cell(1,length(behaviors));
for i = 1:length(behaviors)
    idx = strcmp(event_table.Behavior, behaviors{i}) & strcmp(event_table.BehaviorType, 'START');
    frame_start = event_table.ImageIndex(idx);
    clear tmp


    idx = arr2idx(aggregate(table2array(stroke_events(:,i)),1,fs));
    frame_start = idx(:,1);

    for j = 1:length(frame_start)
        try 
            tmp(:,:,j) = Nresample(:, frame_start(j)-(sz*fs):frame_start(j)+(sz*fs));
            baseline_tmp(:,j) = mean(Nresample(:,frame_start(j)-(sz*fs):frame_start(j)),2);
            response_tmp(:,j) = mean(Nresample(:,frame_start(j):frame_start(j)+(sz*fs)),2);
        catch
            disp('Inside catch')
        end
    end
    
    % Do a paired signrank test for each neuron's baseline vs response to 
    % determine if it is significantly modulated by the behavior
    no_mod_count = 1;
    neg_mod_count = 1;
    pos_mod_count = 1;
    for j = 1:size(baseline_tmp,1)
        [p, h] = signrank(baseline_tmp(j,:), response_tmp(j,:));
%         disp(p)
        if h
            if mean(baseline_tmp(j,:)) > mean(response_tmp(j,:))
                neg_mod{i}(neg_mod_count,:,:) = tmp(j,:,:);
                neg_mod_count = neg_mod_count + 1;
            else
                pos_mod{i}(pos_mod_count,:,:) = tmp(j,:,:);
                pos_mod_count = pos_mod_count + 1;
            end
        else
            no_mod{i}(j,:,:) = tmp(j,:,:);
            no_mod_count = no_mod_count + 1;
        end
    end

    tmp = mean(tmp,3);

    mr = mean(tmp(:,sz*fs:end), 2);
    [~, I] = sort(mr);
%     if i == 5
%         I5 = I;
%     end
    subplot(2, length(behaviors)+2, i)
    imagesc(tmp(I,:))
    hold on
    vline(sz*fs+1, 'k-')
    xticklabels([])
    
    title(behaviors{i})
%     colorbar,
    caxis([-2 2])
    colormap(bluewhitered())

    subplot(2, length(behaviors)+2, i+length(behaviors)+2)
    plot(xt(mean(tmp),fs)-sz, mean(tmp));
    hold on
    axis([-sz sz -0.2 1.5])
    vline(0, 'k')
    xlabel('Time (s)')

end

clear tmp

flrmov = arr2idx(aggregate(flrv >mean(flrv) + std(flrv), 1));
% flrmov = arr2idx(aggregate(lick_events, 3));
% flrmov = arr2idx(flrv >mean(flrv) + std(flrv));
frame_start = flrmov(:,1);
for j = 1:size(flrmov,1)
    try
        tmp(:,:,j) = Nresample(:, frame_start(j)-(sz*fs):frame_start(j)+(sz*fs));
    catch
        disp('in catch 1')
    end
end

tmp = mean(tmp,3);
mr = mean(tmp(:,sz*fs:end), 2);
[~, I] = sort(mr);

subplot(2, length(behaviors)+2, length(behaviors)+1)
imagesc(tmp(I,:))
hold on
vline(sz*fs+1, 'k-')
xticklabels([])
title('FLR move')
%     colorbar,
caxis([-2 2])
colormap(bluewhitered())

subplot(2, length(behaviors)+2, 2*(length(behaviors))+3)
plot(xt(mean(tmp),fs)-sz, mean(tmp));
axis([-sz sz -0.2 1.5])
hold on
vline(0, 'k')

xlabel('Time (s)')



clear tmp

fllmov = arr2idx(aggregate(fllv >mean(fllv) + std(fllv), 1));
frame_start = fllmov(:,1);
for j = 1:size(fllmov,1)
    try
        tmp(:,:,j) = Nresample(:, frame_start(j)-(sz*fs):frame_start(j)+(sz*fs));
    catch
        disp('in catch')
    end
end

tmp = mean(tmp,3);
mr = mean(tmp(:,sz*fs:end), 2);
[~, I] = sort(mr);

subplot(2, length(behaviors)+2, length(behaviors)+2)
imagesc(tmp(I,:))
hold on
vline(sz*fs+1, 'k-')
xticklabels([])
title('FLL move')
%     colorbar,
caxis([-2 2])
colormap(bluewhitered())

subplot(2, length(behaviors)+2, 2*(length(behaviors)+2))
plot(xt(mean(tmp),fs)-sz, mean(tmp));
axis([-sz sz -0.2 1.5])
hold on
vline(0, 'k')

xlabel('Time (s)')


%%

figure
clear tmp
for i = 1:length(behaviors)
    tmp = neg_mod{i};
    
    mr = mean(tmp(:,sz*fs:end,:), [3, 2]);
    [~, I] = sort(mr);

    baseline = mean(tmp(:,1:sz*fs,:), [3 2]);


    subplot(2, length(behaviors), i)
    imagesc(mean(tmp(I,:,:),3) - baseline)
    title(behaviors{i})

    subplot(2, length(behaviors), length(behaviors)+i)
    plot(mean(tmp(I,:,:), [1 3]))
    hold on
    vline(sz*fs, 'k-')
end
caxis([-2 2])
colormap(bluewhitered())

%%

figure,
subplot(1,3,1:2)
imagesc(r), 
caxis([-0.1 0.1])
c=colorbar;
colormap(bluewhitered())
c.Label.String = 'Correlation';
yticklabels(stroke_events.Properties.VariableNames)
xlabel('Neuron Number')
subplot(1,3,3)
xvals = zeros(size(r')) + (1:size(r,1));
% xvals = reshape(xvals, []);
swarmchart(xvals(:), reshape(r', 1, []), 'k', '.')
xticks(1:6)
xticklabels(stroke_events.Properties.VariableNames)



%%

figure,

for j = 1:length(stroke_events.Properties.VariableNames)
    v = stroke_events.Properties.VariableNames{j};
    idx = arr2idx(table2array(stroke_events(:,v)));
    clear tmp
    for i = 1:size(idx,1)
        tmp(:,:,i) = Nresample(:, idx(i,1)-450:idx(i,1)+450);
    end

    [sorted_data, ~] = sort_by_response(mean(tmp,3)');
    subplot(2,ceil(length(stroke_events.Properties.VariableNames)/2), j)
    imagesc(sorted_data'),
    caxis([-2 8])
    colormap(bluewhitered); colorbar
    vline(451, 'k')
    title(v)
end



%%
thresh = 0.01;
figure
for i = 1:size(r,1)
    rcopy = r;
    idx = arr2idx(aggregate(eventsarr(i,:),3,fs));
    clear tmp
%     corr_neurons = r(i,:) > thresh;
    for j = 1:100        
        [rtmp, I] = max(rcopy(i,:));

        for k = 1:size(idx,1)
%                 subplot(2, ceil(size(r,1)/2), i)
%                 plot(xt(1:901, fs)-5, Nresample(I, idx(k,1)-450:idx(k,1)+450), 'color', [0 0 0 0.25]), hold on
%                 title(stroke_events.Properties.VariableNames(i))

            tmp(k,:) = Nresample(I, idx(k,1)-450:idx(k,1)+450);
        end
        subplot(2, ceil(size(r,1)/2), i)
        plot(xt(1:901, fs)-5, mean(tmp), 'color', [0 0 0 0.25]), hold on
        
        rcopy(i,I) = 0;
        vline(0, 'r')
    end
    title(stroke_events.Properties.VariableNames(i))
end

%% plot location of highly correlated neurons onto the brain


I = loadtiff('Y:\nick\behavior\grooming\2p\IDR3_tTA6s\dalsa\AVG_template.tif');
load('Y:\nick\behavior\grooming\2p\IDR3_tTA6s\dalsa\atlas_tform.mat')
load('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat');
test = imwarp(I, tform, 'OutputView', imref2d(size(dorsalMaps.dorsalMapScaled)));

figure, imagesc(test), colormap gray, hold on
axis off;
for p = 1:length(dorsalMaps.edgeOutline)
    plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'w', 'LineWidth', 2);
end
%%
close all
clc
atlas_aligned_fig = fullfile(event_path, '..', 'dalsa', 'atlas_aligned.fig');

uiopen(atlas_aligned_fig,1)
hold on

thresh = 0.15;

for i = 1:size(Nresample,1)
    x0 = double(cstat{i}.med(2));
    y0 = double(cstat{i}.med(1));

    [x,y] = get_neuron_location_in_dorsal_space( x0,y0, ...
        tforms{cstat{i}.use_tform}.roi_to_linear, ...
        tforms{cstat{i}.use_tform}.linear_to_wfield, ...
        tforms{cstat{i}.use_tform}.wfield_to_atlas);

    plot(x, y, 'Color', '#009E73', 'Marker', '.', 'MarkerSize', 5)
end

%%

prev_I = [];
aa=1;
tmp = r(aa,:);

counter = 0;
for i = 1:size(Nresample,1)
    [rr, I] = max(tmp);
    if rr<thresh, continue; end
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
    counter = counter + 1;

end

disp([num2str(counter),' ', stroke_events.Properties.VariableNames{aa},' cells'])

counter = 0;
counter2 = 0;
aa = 2;
tmp = r(aa,:);
for i = 1:size(Nresample,1)
    [rr, I] = max(tmp);
    if rr<thresh, continue; end
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
        counter2 = counter2+1;
    else
        plot(x, y, 'b.', 'MarkerSize', 20)
    end
    counter = counter + 1;
end

disp([num2str(counter),' ', stroke_events.Properties.VariableNames{aa},' cells'])
disp([num2str(counter2),' shared cells'])
% axis tight
% yticklabels(fliplr(nlabel));
% title('Elliptical')






%% plot some example neurons

idx = 7;

t = xt(Nresample, 90);
clc
% num_cells = 20
nlabel = cell(1,10);
figure, hold on%  subplot(1,2,1), hold on
tmp = r(idx,:);
for i = 1:10
    [~, I] = max(tmp);
    plot(t, Nresample(I,:)-10*i, 'k') % "Elliptical" correlated neuron
    tmp(I) = nan;
    nlabel{1,i} = [nloc{I},'-N',num2str(I)];
    
end
axis tight
yticklabels(fliplr(nlabel));
title(stroke_events.Properties.VariableNames{idx})

% vline(t(find(table2array(events(:,1)))), 'r:')

patchplot(t(arr2idx(aggregate(eventsarr(idx,:)', 3))), ylim, 'c', 0.5)
% patchplot(t(arr2idx(aggregate(eventsarr(5,:)',3))), ylim, 'm', 0.5)
xlabel('Time (s)')

%%
close all
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
vline(t(find(table2array(events(:,1)))), 'r:')

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
    rmat = zeros(size(Nresample,1), size(stroke_types,2));
    for i = 1:size(Nresample,1)
        for j = 1:size(stroke_types,2)
            r_tmp = corr(Nresample(i,:)', stroke_types(:,j));
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


%%

all_groom = any(table2array(stroke_events),2);

for i = 1:size(Nresample,1)
    rmov(i) = corr(flrv, Nresample(i,:)');
end

%%


[~, I] = sort(rgroom, 'descend');
figure, 
subplot(1,3,1)
plot(sort(rgroom), 1:length(rgroom))
xlabel('Correlation')
ylabel('Neuron #')
axis tight
subplot(1,3,2:3)
imagesc(imadjust(Nresample(I,:)))
colormap(flipud(colormap('gray')))
patchplot(arr2idx(aggregate(any(table2array(stroke_events(:,[3,5])),2),3)), [0 size(Nresample,1)], 'm', 0.5)
patchplot(arr2idx(aggregate(any(table2array(stroke_events(:,[4,6])),2),3)), [0 size(Nresample,1)], 'c', 0.5)

patchplot(arr2idx(aggregate(any(table2array(stroke_events(:,1)),2),3)), [0 size(Nresample,1)], 'y', 0.5)
patchplot(arr2idx(aggregate(any(table2array(stroke_events(:,2)),2),3)), [0 size(Nresample,1)], 'g', 0.5)


%% unilateral vs bilateral

uni = any(table2array(stroke_events(:,[3,5])),2);
bi = any(table2array(stroke_events(:,[1,2,4,6])),2);

for i = 1:size(Nresample,1)
    runi(i) = corr(uni, Nresample(i,:)');
    rbi(i) = corr(bi, Nresample(i,:)');
end

%%

figure, plot(sort(runi)), hold on, plot(sort(rbi))


%% neural trajectories

[coeff, score, latent] = pca(Nresample);

%%

bb = 3;
% stroke_events(1,:)
b_idx = arr2idx(table2array(stroke_events(:,bb)));

figure, hold on
for i = 1:size(b_idx,1)
    plot(coeff(b_idx(i,1):b_idx(i,2),1), coeff(b_idx(i,1):b_idx(i,2),4), 'r')
end

%%

bb = 1;
% stroke_events(1,:)
b_idx = arr2idx(table2array(stroke_events(:,bb)));

% figure, hold on
for i = 1:size(b_idx,1)
    plot3(coeff(b_idx(i,1):b_idx(i,2),1), coeff(b_idx(i,1):b_idx(i,2),3), coeff(b_idx(i,1):b_idx(i,2),5), 'b')
end


%%
vid_end = find(events.("Video End"));

vel = readmatrix('Y:\nick\behavior\grooming\2p\ECL3_thy1\20240731\ECL3_thy1_20240731_trimDLC_resnet50_mouse_groomingNov27shuffle1_1030000_vel.csv');

flrv = sum(vel(1:vid_end,4:5).^2, 2).^0.5;
fllv = sum(vel(1:vid_end,7:8).^2, 2).^0.5;

%%

[beh, ann] = parse_snippets('Y:\nick\behavior\grooming\2p\RR3_tTA8s\20240802\snippets');
%%

figure, plot(flrv)
hold on, plot(fllv-10)



[coeff, score, latent] = pca(Bmean);
figure
scatter3(score(:,1), score(:,2), score(:,3), 'o')




all_groom = any(table2array(stroke_events),2);
all_groom = aggregate(all_groom, 3);

for i = 1:size(Nresample,1)
    rgroom(i) = corr(all_groom, Nresample(i,:)');
end


%%

[~, I] = sort(rgroom, 'descend');

%%

%% top correlated neurons with grooming
figure, plot(mean(Nresample(I(1:50),:)))

%%

groom_idx = arr2idx(all_groom);


cn = Nresample(I(1:50),:);

%% fit line to each grooming event


fs = 90;
clc
m = [];
leng = [];
angg = [];
for i = 1:size(groom_idx,1)
    for j = 10%:size(cn,1)
        Y = cn(j,groom_idx(i,1):groom_idx(i,2));
    
        X = xt(Y,fs);
        leng = [leng X(end)-X(1)];

        c = polyfit(X, Y, 1);
        m = [m c(1)];
        y_est = polyval(c,X);

        angg = [angg atan((y_est(end)-y_est(1))./(X(end)-X(1)))];

%         figure, plot(X,Y);
%         hold on, plot(X, y_est)

        title(num2str(c(1)))

    
    end

end
%%
[~, aa] = sort(leng);
figure, plot(leng(aa), m(aa))

%% exponential decay (not working well)
fs = 90;
clc
clear tmp
aa = [];
leng = [];
% figure, hold on
for i = 1:size(groom_idx,1)
    for j = 1%:size(cn,1)
        Y = test(groom_idx(i,1):groom_idx(i,2));
    %     plot(tmp-i)
        leng(i) = length(Y);
    %     aa(i) = tmp(end)-tmp(1);
    
        X = xt(Y,fs);
    
    
    
        % Convert X and Y into a table, which is the form fitnlm() likes the input data to be in.
        % Note: it doesn't matter if X and Y are row vectors or column vectors since we use (:) to get them into column vectors for the table.
        tbl = table(X(:), Y(:));
        % Define the model as Y = a * exp(-b*x) + c
        % Note how this "x" of modelfun is related to big X and big Y.
        % x((:, 1) is actually X and x(:, 2) is actually Y - the first and second columns of the table.
        modelfun = @(b,x) b(1) * exp(-b(2)*x(:, 1)) + b(3); 
    
        % Guess values to start with.  Just make your best guess.
        aGuessed = mean(Y); % Arbitrary sample values I picked.
        bGuessed = 0.005;
        cGuessed = Y(1);
        beta0 = [aGuessed, bGuessed, cGuessed]; % Guess values to start with.  Just make your best guess.
        % Now the next line is where the actual model computation is done.
        mdl = fitnlm(tbl, modelfun, beta0);
        % Now the model creation is done and the coefficients have been determined.
        % YAY!!!!
        % Extract the coefficient values from the the model object.
        % The actual coefficients are in the "Estimate" column of the "Coefficients" table that's part of the mode.
        coefficients = mdl.Coefficients{:, 'Estimate'};
        % Create smoothed/regressed data using the model:
        yFitted = coefficients(1) * exp(-coefficients(2)*X) + coefficients(3);
    
        tmp(:,i) = coefficients;
    end

end

%%


[~,aa] = sort(leng);
figure, plot(leng(aa), m(aa))



%%

%% correlate top PCs with behaviors

r = [];
for i = 1:size(bmat,2)
    for j = 1:size(bmat,2)
        r(i,j) = corr(coeff(:,i), bmat(:,j));        
    end
end

figure, imagesc(r), 
c=colorbar;
c.Label.String = 'Correlation'
caxis([-0.2 0.2]);
colormap(bluewhitered())

ylabel('Neural PCs')
yticks(1:size(bmat,2))
xticks(1:size(bmat,2))
xticklabels(regLabels)




