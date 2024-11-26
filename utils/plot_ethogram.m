function plot_ethogram(tab, states, fs)
    figure, hold on
    counter = 1;
    cols = [[0 0.4470 0.7410];
        [0 0.4470 0.7410];
        [0 0.4470 0.7410];
        [0.8500 0.3250 0.0980];
        [0.8500 0.3250 0.0980];
        [0.8500 0.3250 0.0980];
        [0.9290 0.6940 0.1250];
        [0.9290 0.6940 0.1250];
        [0 0.4470 0.7410]];
    t = xt(tab, fs, 1);
    clear myticks
    for i = 1:length(states)
        tmp_idx = find(strcmp(tab.Properties.VariableNames, states(i)));
        if tmp_idx
            if sum(table2array(tab(:,tmp_idx))) > 0
                patchplot(t(arr2idx(table2array(tab(:,tmp_idx)))), [counter counter+1], cols(i,:), 1)
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