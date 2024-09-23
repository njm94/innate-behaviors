addpath('C:\Users\user\Documents\Nick\grooming\utils')
mice = {'Y:\nick\behavior\grooming\1p\ECR2_thy1', ...
    'Y:\nick\behavior\grooming\1p\GER2_ai94', ...
    'Y:\nick\behavior\grooming\1p\HYL3_tTA', ...
    'Y:\nick\behavior\grooming\1p\IBL2_tTA'};

for i = 2%1:length(mice)
    clear mask
    if ~isfile([mice{i}, filesep, 'mask.mat'])
        template = loadtiff([mice{i}, filesep, 'template.tif']);


        mask = draw_roi(template, 2);
%         save([mice{i}, filesep, 'mask.mat'], 'mask')
    end
end