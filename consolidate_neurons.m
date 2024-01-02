function consolidate_neurons(output_file, path_array)
% Consolidate suite2p outputs from separate ROIs into a single mat file
% Input:
%   output_file (path and filename to save consolidated .mat)
%   path_array  (cell array containing paths to suite2p Fall.mat files)


for i=1:length(path_array)
    data = load(path_array{i});
   
    stat{i} = data.stat;
    ops(i) = data.ops;
    F{i} = data.F;
    Fneu{i} = data.Fneu;
    spks{i} = data.spks;
    iscell{i} = data.iscell;
    redcell{i} = data.redcell;
    n_loc{i} = data.n_loc;
end

save(output_file, 'stat', 'ops', 'F*', 'iscell', 'redcell', 'n_loc')
disp('Finished')
end