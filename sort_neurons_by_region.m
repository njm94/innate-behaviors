function [n_sorted, idx, regions] = sort_neurons_by_region(data, nlocs)
idx = [];
n_sorted = [];
regions = {'MO', 'SS', 'RSP', 'VIS'};
for i = 1:length(regions)
    tmp_idx = contains(nlocs, regions{i});
    tmp_neurons = data(tmp_idx, :);
    n_sorted = [n_sorted; tmp_neurons];
    idx = [idx; i+zeros(size(tmp_neurons, 1), 1)];
end

end