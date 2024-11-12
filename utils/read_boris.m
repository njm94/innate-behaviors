function [events, b_idx, t, video_end] = read_boris(boris_tsv, len)
% Reads the annotated BORIS file to extract behaviors
%
% Input:
%       boris_tsv       (Full path to BORIS .tsv events file)
%       len             (Trial recording length (optional). This will use
%                       the last entry in the ImageIndex Column by default.
%                       Specify length in 2 photon trials where end index
%                       is listed in the data_list.py)
%
% Output:
%       events          (Binary table [trial_length x num_behaviors])
%       b_idx           (Cell array of event indices [trial_length x 1])
%       t               (Raw table)
%       video_end       (Ending frame)
%
% Usage: 
%       [events, b_idx, t] = read_boris(boris_tsv);
%       [events, b_idx, t] = read_boris(boris_tsv, len);

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
t = readtable(boris_tsv, "FileType","text",'Delimiter', '\t');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')

if nargin < 2 || isempty(len), len = max(t.ImageIndex); end

% Add 1 to frame index since originally computed with 0-indexing
t.ImageIndex = t.ImageIndex + 1; 

% Initialize the binary event vector
annotations = unique(t.Behavior);
events = zeros(length(annotations), len);
b_idx = cell(1,length(annotations));

% By default set video end to last detected behavior
if any(strcmp(annotations, 'Video End'))
    video_end = t.ImageIndex(strcmp(t.Behavior, 'Video End'));
else
    disp('Video end was not labelled in Boris. Using last behavior as end')
    video_end = t.ImageIndex(end);
end

% Populate event vector
for i = 1:length(annotations)
    behavior_category_index = strcmpi(annotations(i), t.Behavior);
    switch t.BehaviorType{find(behavior_category_index, 1)}
        case 'POINT'
            events(i, t.ImageIndex(behavior_category_index)) = 1;
        case 'START'
            
            behavior_frames = t.ImageIndex(behavior_category_index);
            start_idx = behavior_frames(1:2:end);
            stop_idx = behavior_frames(2:2:end);

            assert(length(start_idx) == length(stop_idx), ...
                'UNEVEN NUMBER OF START/STOP EVENTS...')
            b_idx{i} = [start_idx stop_idx];
            for j = 1:length(start_idx)
                events(i, start_idx(j):stop_idx(j)) = 1;
            end
        otherwise
            disp('ATTENTION: This should not execute, since START always occurs before STOP')
    end
end

% Convert to table
events = array2table(events', 'VariableNames', annotations');

end