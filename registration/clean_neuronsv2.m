clear, clc

[consolidated_file, experiment_path] = uigetfile('*.mat','Select consolidated neuron file.', 'Y:\nick\behavior\grooming\2p');
load([experiment_path, consolidated_file]);
disp('[+] Finished loading consolidated neuron data')

[resampled_data, fs] = resample2maxFs(F, ops);
% group left hemisphere
%%
clc
[X, nloc] = nclean(resampled_data, iscell, n_loc);

%% Clean up stat array

for i = 1:length(stat)
    for j = 1:length(stat{i})
        stat{i}{j}.use_tform = i;
    end
end
cstat = catcell(2, stat)';
good_idx = catcell(1, iscell);
good_idx = logical(good_idx(:,1));
cstat = cstat(good_idx);
% N = [N_left; N_right];
% dF_traces = (N-mean(N,2))./mean(N,2);


% save([experiment_path, 'Fcascade.mat'], 'dF_traces')

%%
X = zscore(X, [], 2);


save([experiment_path, 'Fclean.mat'], 'X', 'nloc', 'fs', 'cstat', 'tforms', '-v7.3')

function [N, nloc] = nclean(data, iscell, n_loc)
% Input:
%       F         (cell array of neurons by time)
%       iscell    (cell array of logicals containing curated neurons)
%
% Output:
%       N         (NxT matrix) 
%       loc       (Nx1) array of regions




iscell = catcell(1, iscell);
nloc = catcell(2, n_loc)';
N = catcell(1, data);

good_idx = iscell(:,1) == 1;

N = N(good_idx, :);
nloc = nloc(good_idx);

end

