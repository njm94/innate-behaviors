function [behaviors, annotations, vPaths] = parse_snippets(snippets_dir)
% Usage:
%   [behaviors, annotations, vPaths] = parse_snippets(snippets_dir)

% annotations = ["Left", "Right", "Elliptical", "LargeLeft", "LargeRight", "LargeBilateral", "Lick"];
annotations = ["elliptical", "largeleft", "largeright", "largebilateral", "left", "right", "lick", 'bilateral'];

snippets = convertCharsToStrings(getAllFiles(snippets_dir, '.tif'));
snippets(contains(snippets, 'discard')) = [];

[behaviors, vPaths] = prune_snippets(snippets, annotations);
end


function [start_idx, end_idx] = get_timestamps_from_str(str)
    expression = '\d+_\d+';
    matchStr = regexp(str,expression,'match');
    
    indices = str2double(split(matchStr, '_'));
    start_idx = indices(1);
    end_idx = indices(2);
end

function [behaviors, vPaths] = prune_snippets(snippets, annotations)
    behaviors = cell(1, length(annotations));
    vPaths = cell(1, length(annotations));
    for i = 1:length(annotations)
        if ~any(contains(snippets, annotations(i)))
            continue
        else
            index = contains(snippets, annotations(i));
            vPaths{i} = snippets(index);
            for j = 1:length(vPaths{i})
                [behaviors{i}(j,1), behaviors{i}(j,2)] = get_timestamps_from_str(vPaths{i}{j});
            end
            snippets(index) = [];
        end
    end

    index = contains(snippets, "discard");
    snippets(index) = [];

    if ~isempty(snippets)
        disp(['Improperly labeled behavior annotation was detected: ', snippets{:}])
    end
end