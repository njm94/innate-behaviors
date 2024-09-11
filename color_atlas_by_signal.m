% Fill in atlas with mean dFF values

clear, clc

experiment_root_dir = 'Y:\nick\behavior\grooming\2p';
mouse_id = 'ECL3_thy1';

dc = dir([experiment_root_dir, filesep, mouse_id]);
dc = {dc.name};
ignore = contains(dc, {'dalsa', 'outputs', '.', 'ignore'});

dc = dc(~ignore);
all_possible_behaviors = {'Elliptical', 'Elliptical Asymmetric', 'Left', ...
    'Right', 'Left Asymmetric', 'Right Asymmetric', 'Lick Unkown', ...
    'Large Bilateral'};

regions = cell(1,length(all_possible_behaviors));
neuron_signals = cell(1,length(all_possible_behaviors));
for i = 1:length(dc)
    
    if isempty(getAllFiles([experiment_root_dir, filesep, mouse_id, filesep, dc{i}], '.tsv'))
        disp(['Skipping ', dc{i}, '... No TSV found.'])
        continue
    else
        boris_file = getAllFiles([experiment_root_dir, filesep, mouse_id, filesep, dc{i}], '.tsv');
        [events, b_idx, ~] = read_boris([experiment_root_dir, filesep, mouse_id, filesep, dc{i}, filesep, boris_file]);
        
        if ~isfile([experiment_root_dir, filesep, mouse_id, filesep, dc{i}, filesep, 'Nresample.mat'])
            load([experiment_root_dir, filesep, mouse_id, filesep, dc{i}, filesep, 'Fclean.mat']);
            
            disp("Resampling neural data to match behavior")
            try Nresample = resample(N, size(events,1), size(N, 2), 'Dimension', 2);
            catch
                disp('Data too big. Splitting into halves and resampling each half separately')
                Nresample = resamplee(N', size(events,1), size(N,2))';
            end
        
            save([experiment_root_dir, filesep, mouse_id, filesep, dc{i}, filesep, 'Nresample.mat'], ...
                'Nresample', 'cstat', 'nloc', 'tforms')
        else
            disp('Loading resampled neuron data')
            load([experiment_root_dir, filesep, mouse_id, filesep, dc{i}, filesep, 'Nresample.mat'])
        end

        % load DLC tracks
        vel_file = getAllFiles([experiment_root_dir, filesep, mouse_id, filesep, dc{i}], '_vel.csv');
        vel = readmatrix([experiment_root_dir, filesep, mouse_id, filesep, dc{i}, filesep, vel_file]);
        vid_end = find(events.("Video End"));        
        vel = vel(1:vid_end,:);
        flrv = sum(vel(:,4:5).^2, 2).^0.5;
        fllv = sum(vel(:,7:8).^2, 2).^0.5;       
        flrthresh = flrv>mean(flrv) + std(flrv);
        fllthresh = fllv>mean(fllv) + std(fllv);

        for j = 1:length(all_possible_behaviors)
%             if ~strcmpi(events.Properties.VariableNames, all_possible_behaviors{j})
%                 continue
%             else
                behavior_idx = strcmpi(events.Properties.VariableNames, all_possible_behaviors{j});
                Bmean = mean(Nresample(:,  logical(table2array(events(:,behavior_idx)))),2);

                regions_from_current_session = unique(nloc);
                for k = 1:length(regions_from_current_session)
                    if any(strcmpi(regions{j}, regions_from_current_session{k}))
                        idx = find(strcmpi(regions{j}, regions_from_current_session{k}));
                        neuron_signals{j}{idx} = [neuron_signals{j}{idx}; Bmean(contains(nloc, regions_from_current_session{k}))];
                    else
                        regions{j} = [regions{j}; regions_from_current_session(k)];
                        neuron_signals{j} = [neuron_signals{j}; {Bmean(strcmpi(nloc, regions_from_current_session{k}))}];
                    end

                end
%             end
        end




    end
end


%%


load('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat');
load('C:\Users\user\Documents\Nick\grooming\utils\atlas.mat')
         
figure
for i = 1:length(all_possible_behaviors)
    subplot(2,4,i)
    

end


%%
