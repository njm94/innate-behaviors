% Consolidate suite2p outputs from separate ROIs into a single mat file
% Input:
%   output_file (path and filename to save consolidated .mat)
%   path_array  (cell array containing paths to suite2p Fall.mat files)


clc, clear

startpath = 'Y:\nick\2p';
working_dir = uigetdir(startpath, 'Select experiment directory.');

arm_array = get_folders(working_dir);

for i = 1:length(arm_array)
    roi_array{i} = get_folders(arm_array{i});
    for j = 1:length(roi_array)
        plane_array{j} = get_folders(roi_array{i});
    end
end


roi_array = catcell(2, roi_array);
template_path = get_template_path(working_dir);

for i=1:length(roi_array)
    for j = 1:length(plane_array)
    data = load([roi_array{i}, '\bin2x2x1\suite2p\plane0\Fall.mat']);
   
    stat{i} = data.stat;
    ops(i) = data.ops;
    F{i} = data.F;
    Fneu{i} = data.Fneu;
    spks{i} = data.spks;
    iscell{i} = data.iscell;
    redcell{i} = data.redcell;
    n_loc{i} = data.n_loc;


    roi_to_linear = load([roi_array{i}, '\bin2x2x1\tform.mat']);


    % if there are separate linear scans for each arm load them
    % otherwise, load linear scan from base working directory (date)
    pathparts = strsplit(roi_array{i}, filesep);
    arm_path = join(pathparts(1:end-1), filesep);
    try
        linear_to_wfield = load([arm_path{:}, filesep, 'linearscan_tform.mat']);
    catch
        linear_to_wfield = load([working_dir, filesep, 'linearscan_tform.mat']);
    end
    
    wfield_to_atlas = load([template_path, filesep, 'atlas_tform.mat']);
    tforms{i}.roi_to_linear = roi_to_linear;
    tforms{i}.linear_to_wfield = linear_to_wfield;
    tforms{i}.wfield_to_atlas = wfield_to_atlas;
    end
    
end

save([working_dir, filesep, 'Fall.mat'], 'stat', 'ops', 'F*', 'iscell', 'redcell', 'n_loc', 'tforms', 'spks')
disp('Finished')


function path_array = get_folders(working_dir)

clist = dir(working_dir);
flist = [clist.isdir];

count = 1;
for i = 1:length(flist) 
    if flist(i) && ~contains(clist(i).name, '.')
        path_array{count} = [clist(i).folder, filesep, clist(i).name];
        count = count + 1;
    end
end


end


function template_path = get_template_path(working_dir)

folders = strsplit(working_dir, '\');
for i = length(folders):-1:1
    tmp = join(folders(1:i), '\');
    tmp = get_folders(tmp{:});
    if any(contains(tmp, 'template'))
        template_path = tmp{contains(tmp, 'template')};
        break
    end

end

end