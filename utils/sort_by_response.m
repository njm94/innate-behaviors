function [sorted_data, I] = sort_by_response(data)
% data is Time x Neuron x Trial
%

avg_data = mean(data, 3); % average across trials
stim_time = round(size(data,1)/2);

[~, tmpI] = max(avg_data(stim_time:end, :));
[~, I] = sort(tmpI);
sorted_data = data(:, I, :);