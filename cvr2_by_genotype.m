


fs = 90;
fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
data_list = textscan(fileID, formatSpec);


for j = 6:length(data_list{1})
    data_dir = data_list{1}{j};
    disp(['Starting ' data_dir])
    fPath = [data_dir filesep 'ridge_outputs_ipsi_contra_bilatInc_forelimbMovExc_audio_drop'];
    try
        load([fPath, filesep, 'cvFull.mat'])
    catch
        continue
    end

    disp('Loading brain data...')
    brain_file = [data_dir, filesep, get_file_with_str(data_dir, 'cam0_svd')];
    load(brain_file)
    Ubrain = U;
    Vbrain = s*V(:,1:end-5);
    [b, a] = butter(2, 0.01/(fs/2), 'high');
    Vbrain = filtfilt(b, a, Vbrain')';

    if size(Vbrain,2) < size(Vfull,2)
        Vfull = Vfull(:,1:end-5);
    end

    fullMat = modelCorr(Vbrain,Vfull,Ubrain) .^2;
    cvr2_avg(j) = mean(fullMat);

end


%%
thy1 = cvr2_avg(1:5)';
ai94 = cvr2_avg(8:12)';
camk = cvr2_avg([14:18,20:24])';

figure,
uboxplot(thy1, ai94, camk)
hold on, plot(ones(size(thy1))+0.3, thy1, 'ko')
plot(ones(size(ai94))+1.3, ai94, 'ko')
plot(ones(size(camk))+2.3, camk, 'ko')
ylabel('Explained Variance (Mean over cortex)')
xticklabels({'Thy1', 'Ai94', 'CamKII'})


