clc, clear
fileID = fopen('expt3_datalist.txt','r');
addpath('C:\Users\user\Documents\Nick\grooming\utils')

formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
data_list = data_list{1};

hyl3_idx = 52:56;
ecr2_idx = 57:59;

%%
clc
for i = [hyl3_idx(2)]
    fpath = data_list{i};
    [events, b_idx, b_tab] = read_boris([fpath, filesep, getAllFiles(fpath, 'events.tsv')]);
    load([fpath, filesep, getAllFiles(fpath, 'cam0_svd.mat')])
end

%%

ll = zeros(1,size(V,2));
for i = 1:size(b_idx{1},1)
ll(b_idx{1}(i,1):b_idx{1}(i,2)) = 1;
end

%%

fs = 30;
[b, a] = butter(2, 0.01/(fs/2), 'high');

% test = filtfilt(b, a, V(1,:));

t = xt(V, fs, 2);
figure, plot(t, -V(1,:)), hold on, patchplot(t(b_idx{1}), [-8e-2 8e-2], 'm')





