function path = fix_path(path)
% This function corrects the path if working on the Linux computer or
% Windows

if ~isempty(dir(path)), return; end

if isunix
    path = strrep(path, 'Y:', '/media/user/teamshare');
    path = strrep(path, '\', '/');
else
    path = strrep(path, '/media/user/teamshare', 'Y:');
    path = strrep(path, '/', '\');
end