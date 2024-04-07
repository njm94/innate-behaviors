%% align mice to common template
clear, clc


data_root = 'Y:\nick\behavior\grooming\1p';
mice = {'ECR2_thy1', 'GER2_ai94', 'HYL3_tTA', 'IBL2_tTA'};

template = zeros(128, 128, numel(mice));

for i = 1:length(mice)
    template(:,:,i) = loadtiff([data_root, filesep, mice{i}, filesep, 'template.tif']);
end

%% select bregma for all mice 

for i = 1:size(template, 3)
    figure, imagesc(template(:,:,i)), colormap gray, grid on
    [x(i), y(i)] = ginput(1); 
    close(gcf)
end

%% using ecr2 as template, align all the others to that one


for i = 1:4

    tx = x(1)-x(i);
    ty = y(1)-y(i);
    tform = transltform2d(tx,ty);


    test(:,:,i) = imwarp(template(:,:,i), tform, 'OutputView', imref2d(size(template(:,:,i))));
end
