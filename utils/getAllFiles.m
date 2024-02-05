function fileList = getAllFiles(dirName, ext)
% return list of files from specified directory 
%
% input:
%   dirName (directory name, string)
%   ext     (string, eg .tif, .mat, .csv, etc)
%
% output:
%   fileList (list of file names, cell array of strings)

if nargin < 2 || isempty(ext), ext = '.'; end


dirData = dir(dirName); 

% exclude index for directories
dirIndex = [dirData.isdir]; 
fileList = {dirData(~dirIndex).name}'; 

% restrict fileList to specified filename extension
fileList = fileList( contains(fileList, ext) );

if length(fileList) == 1, fileList = cell2mat(fileList); end

end