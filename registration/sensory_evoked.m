function [dFF, fluo, ref] = sensory_evoked(is_strobing, do_filter, save_raw)
% Trial-based sensory-evoked analysis for mesoscale imaging.
% This assumes that one stimulus occurs half-way through the recording and
% that each file has the same number of frames.
% 
%

if nargin<1 || isempty(is_strobing); is_strobing = true; end
if nargin<2 || isempty(do_filter); do_filter = false; end
if nargin<3 || isempty(save_raw); save_raw = true; end

warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

% select files
[file,path] = uigetfile('*.tif',...
   'Select One or More Files', ...
   'MultiSelect', 'on');
dFF = [];
dFF1 = [];
dFF2 = [];

for i=1:length(file)
    if iscell(file) == 1
        filename = file{i};
    else
        filename = file;
        if i<length(filename)
            continue
        end
    end
    I=loadtiff([path, filename]);
    
    if is_strobing
        % de-interleave
        I1 = I(:,:,1:2:end);
        I2 = I(:,:,2:2:end);
        
        % calculate dFF
        % baseline period starts at 2 in case of dim first frame
        baseline_end = round(size(I1,3)/2);
        tmp1 = calc_dFF(I1, [2, baseline_end]); 
        tmp2 = calc_dFF(I2, [2, baseline_end]);  
        
        % determine fluorescence and reflectance by comparing std of the
        % time-series (assume SD(fluo) > SD(ref) )
        SD1 = mean(mean(std(tmp1,[],3)));
        SD2 = mean(mean(std(tmp2,[],3)));
        if SD1>SD2
            fluo = tmp1;
            ref = tmp2;
        else
            fluo = tmp2;
            ref = tmp1;
        end
        
        % correct for hemodynamics by subtracting reflectance signal
        dFF_corrected = fluo-ref;
    
    else
        % calculate dFF
        % baseline period starts at 2 in case of dim first frame
        dFF_corrected = calc_dFF(I, [2, round(size(I,3)/2)]);
    end
    
    % filter
    if do_filter
        dFF_corrected = image_filter(dFF_corrected);
    end
    
    % concatenate dFF 
    dFF = cat(4, dFF, dFF_corrected);
    
    if is_strobing && nargout > 1
        dFF1 = cat(4, dFF1, fluo);
        dFF2 = cat(4, dFF2, ref);
    end
    
end

% calculate the average dFF
dFF_avg = mean(dFF, 4)*100;

save([path, filesep, 'dFF.mat'], 'dFF')


% save as raw
if save_raw
    fileID = fopen([path, filesep, filename(1:end-3), 'raw'],'w', 'l');
    fwrite(fileID, single(dFF_avg),'single');
    fclose(fileID);
end


end
% 
% function parse_filename(filename)
% exp_date = regexp(filename, '\d{8}', 'match');
% fps = regexp(filename, '\d\d', 'match');
% end

function I = loadtiff(tiff_file)

info = imfinfo(tiff_file);
tampon = imread(tiff_file,'Index',1);
F = length(info);
I = zeros(size(tampon,1),size(tampon,2),F,'uint16');
I(:,:,1) = tampon(:,:,1);
tic
wait_bar = waitbar(0,['Loading ',tiff_file]);
ind = 0;
for i = 2:F
    if ind == 0, waitbar(i/F, wait_bar); end
    ind = ind + 1; if ind == 100, ind = 0; end
    tampon = imread(tiff_file,'Index',i,'Info',info);
    I(:,:,i) = tampon(:,:,1);
end
close(wait_bar);
temps = num2str(round(10*toc)/10);
disp([tiff_file ' open in ' num2str(temps) 's'])
end


function dFF = calc_dFF(F,window)
% Calculate pixel-wise fractional change in fluorescence signal (delta F/F)
% over 3rd dimension
% 
% Inputs:
%   F          (required, image data - dimensions: Height x Width x Time)
%   window     (optional, two element vector -[start end] of baseline period)
%
% Outputs:
%   dFF        (pixel-wise dFF)
% 
% Usage:  dFF = calc_dFF(F, window);

if nargin<2, window(1) = 1; window(2) = size(F,3); end

F = double(F);
F0 = mean(F(:,:,window(1):window(2),:),3);
F0 = repmat(F0,1,1,size(F,3)); 

dFF = (F-F0)./F0;

end


function filt = image_filter(I)
% assume time on 3rd dimension


fs = 150;




[b,a] = butter(2,[31 74]/(fs/2), 'bandpass');
% [b,a] = butter(2,[4 16]/(fs/2), 'bandpass');
%[b,a] = cheby1(3, 5, [4 6]/(fs/2), 'stop');
I = permute(I, [3, 1, 2]);
mp = round(size(I,1)/2);

I = [flipud(I(1:mp,:,:)); I; flipud(I(mp:end,:,:))];
filt = filtfilt(b,a,I);
filt = permute(filt(mp+1:end-mp,:,:), [2, 3, 1]);
end