%% alignment over days
 % Align all trials by concatenating all spatial SVD components from each
 % trial after registration to template, then recompute SVD on master set
 % and retain top 1000 components


clear, close all
addpath('C:\Users\user\Documents\Nick\ridgeModel');
addpath('C:\Users\user\Documents\Nick\ridgeModel\widefield')
addpath('C:\Users\user\Documents\Nick\ridgeModel\smallStuff') 

%%
clc
fileID = fopen('expt1_datalist.txt','r');
formatSpec = '%s';
data_list = textscan(fileID, formatSpec);
current_mouse = '';
for j = 1:length(data_list{1})
    data_dir = data_list{1}{j};
    disp(['Starting ' data_dir])

    [mouse_root_dir, exp_date, ~] = fileparts(data_dir);
    [~, mouse_id, ~] = fileparts(mouse_root_dir);
    [expt_root_dir, ~, ~] = fileparts(mouse_root_dir);
    
%     master_SVD_file = [mouse_root_dir filesep 'masterSVD.mat'];
%     if ~isfile(master_SVD_file)
    
    tform_file = [data_dir filesep 'tform.mat'];
%     if ~isfile(tform_file)
        if ~strcmp(mouse_id, current_mouse)
            if ~isempty(current_mouse)
                % if all the trials have completed for the last mouse,
                % re-compute SVD on the set of concactenated spatial
                % components to obtain the final Umaster
                disp('Running svd on basis set')
                [Umaster, ~, ~] = svd(Umaster);
                Umaster = Umaster(:,1:1000);
                save([expt_root_dir filesep current_mouse filesep 'Umaster.mat'], 'Umaster')
                
            end
            current_mouse = mouse_id;
            mouse_ref_image = loadtiff([mouse_root_dir filesep 'template.tif'], false);
            expt_ref_image = loadtiff([data_dir, filesep, 'cam0_singleFrame.tif'], false);
            Umaster = [];
        end

        tformEstimate = imregcorr(expt_ref_image, mouse_ref_image);
        disp('Saving transformation')
        save(tform_file, 'tformEstimate');
%     else
%         disp('Loading transformation')
%         load(tform_file)
%     end

    try
        brain_file = [data_dir, filesep, 'cam0_svd'];
    catch
        disp('Missing one or more of the required files. Skipping...')
        continue
    end

    % load brain svd components
    disp('Loading brain data...')
    load(brain_file, 'U');
    % Image data needs to be transposed to deal with difference in MATLAB
    % column-major vs Python row-major alignment
    Ubrain = permute(reshape(U, 128, 128, []), [2 1 3]);
    Ubrain = evaluate_tform(Ubrain, tformEstimate); % apply registration
    Ubrain = reshape(Ubrain, 128*128, []);
    Umaster = cat(2, Umaster, Ubrain);
    
%     if j == 5
%         break
%     end

%     fPath = [data_dir filesep 'ridge_outputs_ipsi_contra_bilatInc_forelimbMovExc_audio_drop' filesep];
%     if ~isdir(fPath), mkdir(fPath); end
% %     if isfile([fPath, 'summary.fig'])
% %         disp([data_dir, ' already processed. Skipping...'])
% %         continue
% %     end
% 
%     
%     exp_date = strfind(data_dir, filesep);
%     exp_date = data_dir(exp_date(end)+1:end);
%     fs = 90;
% 



end