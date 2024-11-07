mice = {'Y:\nick\behavior\grooming\1p\ECR2_thy1', ...
    'Y:\nick\behavior\grooming\1p\GER2_ai94', ...
    'Y:\nick\behavior\grooming\1p\HYL3_tTA'};

for i = 1:length(mice)
    cv_file = getAllFiles([mice{i}, filesep, 'outputs'], 'cvFull');
    if size(cv_file,1)>1, cv_file = cv_file{end}; end
    load([mice{i}, filesep, 'outputs', filesep, cv_file])
    load([mice{i}, filesep, 'Umaster.mat'])
    load([mice{i}, filesep, 'mask.mat'])

    cBetaRight = check_beta('Timer', fullLabels, fullIdx, Umaster, fullBeta{1}, Vfull, [], false);

    timer_kernel(:,:,i) = cBetaRight.*mask;
end

%%
timer_kernel(timer_kernel==0) = nan;
figure, hold on
for i = 1:size(timer_kernel, 3)
    tmp = timer_kernel(:,:,i);
    histogram(tmp(:), 50)
end

vline(0, 'k-')


xlabel('Pixel Weights')
ylabel('Count')