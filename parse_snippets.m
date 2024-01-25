function [behaviors, annotations] = parse_snippets(snippets_dir)

annotations = ["lick", "elliptical", "largeleft", "largeright", "largebilateral", "right", "left"];

snippets = getAllFiles(snippets_dir, '.tif');
snippets(contains(snippets, '_discard')) = [];

behaviors = prune_snippets(snippets, annotations);
end


function [start_idx, end_idx] = get_timestamps_from_str(str)
    expression = '\d+_\d+';
    matchStr = regexp(str,expression,'match');
    
    indices = str2double(split(matchStr, '_'));
    start_idx = indices(1);
    end_idx = indices(2);
end

function behaviors = prune_snippets(snippets, annotations)
    behaviors = cell(1, length(annotations));
    for i = 1:length(annotations)
        index = contains(snippets, annotations(i));
        behavior_i = snippets(index);
        for j = 1:length(behavior_i)
            [behaviors{i}(j,1), behaviors{i}(j,2)] = get_timestamps_from_str(behavior_i{j});
        end
        snippets(index) = [];
    end
    if ~isempty(snippets)
        error('Improperly labeled behavior annotation was detected:')
    end
end