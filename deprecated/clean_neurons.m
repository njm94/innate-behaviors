clear, clc

% CHANGE THESE FOR EACH IMAGING SESSION
% Index depends on roi placement
%
% arm1 is always read first in the consolidation script, so indexing starts
% from arm 1 [roi1:N] then goes to arm2 [roi1:N]
left_hem = [];
right_hem = [1, 2, 3];

%% Don't change from here on out

[consolidated_file, experiment_path] = uigetfile('*.mat','Select consolidated neuron file.', 'Y:\nick\behavior\grooming\2p');
load([experiment_path, consolidated_file]);
disp('[+] Finished loading consolidated neuron data')

[resampled_data, fs] = resample2maxFs(F, ops);
% group left hemisphere
%%
clc

[N_left, nloc_left] = nclean(resampled_data, iscell, n_loc, left_hem, false);

[N_right, nloc_right] = nclean(resampled_data, iscell, n_loc, right_hem, true);

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
N = zscore([N_left; N_right], [], 2);
nloc = [nloc_left; nloc_right];

save([experiment_path, 'Fclean.mat'], 'N', 'nloc', 'fs', 'cstat', 'tforms', '-v7.3')

function [N, nloc] = nclean(data, iscell, n_loc, hem_idx, is_right)
% Input:
%       F         (cell array of neurons by time)
%       iscell    (cell array of logicals containing curated neurons)
%       hem_idx   (index of regions to group together)
%       is_right  (bool: do these regions correspond to right hemisphere?)
%
% Output:
%       N         (NxT matrix) 
%       loc       (Nx1) array of regions


empty_idx = find(cellfun(@isempty, iscell));

if isempty(hem_idx)
    N = [];
    nloc = [];
    return
end

iscell = catcell(1, iscell(hem_idx));
nloc = catcell(2, n_loc(hem_idx))';

% In rare cases, no neurons are found in roi. This is to handle those cases
% keeping track of which rois are empty, then decrementing the hemisphere
% index appropriately, since the resampled neuron data is only of length N
% where N is the number of ROIs that do contain neurons

for i = 1:length(hem_idx)
    hem_idx(i) = hem_idx(i) - sum(hem_idx(i)>=empty_idx);
end

hem_idx = hem_idx(hem_idx>0);
disp(hem_idx)
N = catcell(1, data(hem_idx));


good_idx = iscell(:,1) == 1;

N = N(good_idx, :);
nloc = nloc(good_idx);

if is_right
    nloc = cellfun(@(x) strrep(x,'_L','_R'), nloc, 'UniformOutput', false);
end

end

