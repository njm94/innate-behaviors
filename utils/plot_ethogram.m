function plot_ethogram(tab, states, fs)
    % figure, hold on
    counter = 1;
    cols = [[0 0.4470 0.7410]; % Start
        [0 0.4470 0.7410]; % Right
        [0 0.4470 0.7410]; % Left
        [0.8500 0.3250 0.0980]; % Ellip
        [0.8500 0.3250 0.0980]; % R Asymm
        [0.8500 0.3250 0.0980]; % L Asymm
        [0.9290 0.6940 0.1250]; % Ellip R
        [0.9290 0.6940 0.1250]; % Ellip L
        [0 0.4470 0.7410] % Stop
        [1 0 1]; % Lick
        [0 0 0]]; % Drop
    t = xt(tab, fs, 1);
    clear myticks
    for i = 1:length(states)
        if i == 11 % This is the drop. 
            % use contains isntead of strcmp to consolidate dropright
            % dropleft, dropcenter
            tmp_idx = find(contains(tab.Properties.VariableNames, states(i)));
        else
            tmp_idx = find(strcmp(tab.Properties.VariableNames, states(i)));
        end

        if tmp_idx
            if sum(table2array(tab(:,tmp_idx))) > 0
                if i == 11
                    % Drop is a point event, so can't use patchplot
                    % Make drop of length 3 by adding a 1 immediately
                    % before the events
                    tmp_data = table2array(tab(:,tmp_idx));
                    data2patch = t(arr2idx(tmp_data + circshift(tmp_data,-1)+ circshift(tmp_data,-2)));
                else
                    data2patch = t(arr2idx(table2array(tab(:,tmp_idx))));
                end
                patchplot(data2patch, [counter counter+1], cols(i,:), 1)
                myticks(counter) = counter+0.5;
                stateidx(counter) = tmp_idx;
                counter = counter + 1;
            end     
        end
    end
    yticks(myticks)
    yticklabels(tab.Properties.VariableNames(stateidx))
    xlabel('Time (s)')
end