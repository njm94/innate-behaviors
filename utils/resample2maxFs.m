function [resampled_data, new_fs] = resample2maxFs(F, ops)
    [len, idx] = max(cellfun('size', F, 2));
    resampled_data = cell(length(F), 1);
    for i = 1:length(F)
        resampled_data{i} = resample(double(F{i}), len, size(F{i}, 2), 'Dimension', 2);
    end
    
    new_fs = ops(idx).fs;
end