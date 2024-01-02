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