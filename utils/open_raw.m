function I = open_raw(filename, precision, use_blue, resolution)
if nargin < 2 || isempty(precision), precision = 'single'; end
if nargin < 3 || isempty(use_blue); use_blue = false; end
if nargin < 4 || isempty(resolution), resolution = [256, 256]; end

file = dir(filename);
fid = fopen([file.folder,'\',file.name], 'r');

if fid == -1, error(['Can not open ', file.name]); end


disp(['Opening ', filename,' ...'])
I = fread(fid, file.bytes, precision);
fclose(fid);

if strcmp(precision, 'uint8')
    I = uint8(reshape(I, 3, resolution(1), resolution(2), []));
    if ~use_blue
        I([1, 3], :, :, :) = [];
        I = squeeze(I);
    else
        I(1, :, :, :) = [];
        I = permute(I, [2, 3, 4, 1]);
    end
else
    I = reshape(I,resolution(1),resolution(2),[]);
end



% I = permute(I,[2 3 4 1]);

