function [seeds, labels] = get_seeds()

if ~isunix
    load('Y:\nick\2p\code\utils\allen_map\allenDorsalMap.mat');
else
    load('/home/user/Documents/grooming/utils/allen_map/allenDorsalMap.mat')
end
midline = round(size(dorsalMaps.dorsalMapScaled, 2)/2);


figure
hold on
for p = 1:length(dorsalMaps.edgeOutline)
plot(dorsalMaps.edgeOutline{p}(:, 2), dorsalMaps.edgeOutline{p}(:, 1), 'k');
end
set(gca, 'YDir', 'reverse');
vline(midline)

% [x, y] = ginput(5);

seeds = [223 153; % MOS_1
    259 202; % MOS_2
    268 256; % MOS_3
    164 204; % MOP_1
    203 242; % MOP_2
    153 260; % SSP-ul
    197 298; % SSP-ll
    121 237; % SSP-m
    109 277; % SSP-n
    114 337; % SSP-bfd
    262 356; % RSP_1
    251 407; % RSP_2
    169 377; % PTLp
    205 394; % VIS-am
    153 428; % VIS-p
    ];

labels = {'MOS_1', 'MOS_2', 'MOS_3', 'MOP_1', 'MOP_2', 'SSP-ul', ...
    'SSP-ll', 'SSP-m', 'SSP-n', 'SSP-bfd', 'RSP_1', 'RSP_2', 'PTLp', ...
    'VIS-am', 'VIS-p'};

labelsR = cellfun(@(x) cat(2, x, '-R'), labels, 'UniformOutput', false);
labels = [cellfun(@(x) cat(2, x, '-L'), labels, 'UniformOutput', false), labelsR]';
seedsR = seeds;
seedsR(:,1) = midline+midline-seedsR(:,1);

seeds = [seeds; seedsR];

for p = 1:size(seeds,1)
    plot(seeds(p,1), seeds(p,2), 'r*')
end


end