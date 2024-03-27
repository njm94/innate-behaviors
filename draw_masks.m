mice = {'Y:\nick\behavior\grooming\1p\ECR2_thy1', ...
    'Y:\nick\behavior\grooming\1p\GER2_ai94', ...
    'Y:\nick\behavior\grooming\1p\HYL3_tTA', ...
    'Y:\nick\behavior\grooming\1p\IBL2_tTA'};

for i = 1:length(mice)
    clear mask
    if ~isfile([mice{i}, filesep, 'mask.mat'])
        template = loadtiff([mice{i}, filesep, 'template.tif']);
        if contains(mice{i}, 'ai94') % retain olfactory bulbs for ai94
            num_rois = 4;
        else
            num_rois = 2;
        end

        mask = draw_roi(template, num_rois);
        save([mice{i}, filesep, 'mask.mat'], 'mask')
    end
end