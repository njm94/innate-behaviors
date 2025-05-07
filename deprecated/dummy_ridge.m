% Parameters
fs = 100;  % Sampling frequency (Hz)
T = 1000;      % Duration (seconds)
t = 0:1/fs:T-1/fs;  % Time vector
f = 1;     % Frequency of sine wave (Hz)
A = 3;      % Amplitude of sine wave
noise_level = 0.5;  % Standard deviation of noise

% Generate a 1-s sine wave
tt = 0:1/fs:1-1/fs;
sine_wave = A * sin(2 * pi * f * tt);

% Generate a 1-s ramp
ramp = (1:fs)/20;


rng(108);

% Generate a random step sequence
steps = repelem(randperm(5),20);
steps = 3*zscore(steps);

decay = 3*exp(-0.05*(0:99));



noisy_signal = randn(1,T);
figure, 

subplot(5,4,1:4), plot(noisy_signal,'k'), title('Noisy signal')
subplot(5,4,5), plot(sine_wave, 'b'), title('Sine wave')
subplot(5,4,6), plot(ramp,'r'), title('Ramp')
subplot(5,4,7), plot(steps, 'g'), title('Steps')
subplot(5,4,8), plot(decay, 'c'), title('Decay')

sin_idx = randperm(T, 3);
ramp_idx = randperm(T, 3);
step_idx = randperm(T, 3);
decay_idx = randperm(T,3);


% Add Gaussian noise
noisy_signal = addin(noisy_signal, sine_wave, sin_idx);
noisy_signal = addin(noisy_signal, ramp, ramp_idx);
noisy_signal = addin(noisy_signal, steps, step_idx);
noisy_signal = addin(noisy_signal, decay, decay_idx);

subplot(5,4,9:12), plot(noisy_signal, 'k')
hold on, vline(sin_idx, 'b-')
vline(ramp_idx, 'r-')
vline(step_idx, 'g-')
vline(decay_idx, 'c-')

%%
  % ridge here
disp('Building design matrix')
opts.frameRate = fs;
opts.sPostTime=round(fs*1)-1;
opts.framesPerTrial = T; % nr. of frames per trial
opts.folds = 10; %nr of folds for cross-validation

regressor_mat = zeros(T,4);
regressor_mat(ramp_idx,1) = 1;
regressor_mat(sin_idx,2) = 1;
regressor_mat(step_idx,3) = 1;
regressor_mat(decay_idx,4) = 1;
% Full-Trial events:    new_trial
% Peri-Stimulus events: 
[dMat, regIdx] = makeDesignMatrix(regressor_mat, [2, 2, 2, 2], opts);
regLabels = {'Ramp', 'Sine', 'Steps', 'Decay'}; %some movement variables
%             [dMat, regIdx] = makeDesignMatrix(regressor_mat, [3, 3, 1], opts);
%             regLabels = {'Movment', 'Lick', 'Trial'}; %some movement variables
fullR = dMat;

disp('Running ridge regression with 10-fold cross-validation')
[Vfull, fullBeta, ~, fullIdx, fullRidge, fullLabels] = crossValModel(fullR, noisy_signal, regLabels, regIdx, regLabels, opts.folds);
%             save([fPath 'cvFull.mat'], 'Vfull', 'fullBeta', 'fullR', 'fullIdx', 'fullRidge', 'fullLabels', '-v7.3'); %save some results


%%

sine_avg = get_avg(noisy_signal, sin_idx, 100);
ramp_avg = get_avg(noisy_signal, ramp_idx, 100);
steps_avg = get_avg(noisy_signal, step_idx, 100);
decay_avg = get_avg(noisy_signal, decay_idx, 100);


betas = mean(catcell(2, fullBeta),2);


subplot(5,4,13:16), plot(noisy_signal, 'k', 'LineWidth', 1), hold on, plot(Vfull), title('Model reconstructions')
subplot(5,4,17), plot(sine_wave, 'k', 'LineWidth', 1), hold on, plot(betas(101:200)), plot(sine_avg), 
subplot(5,4,18), plot(ramp, 'k', 'LineWidth', 1), hold on, plot(betas(1:100)), plot(ramp_avg),
subplot(5,4,19), plot(steps, 'k', 'LineWidth', 1), hold on, plot(betas(201:300)), plot(steps_avg), 
subplot(5,4,20), plot(decay, 'k', 'LineWidth', 1), hold on, plot(betas(301:400)), plot(decay_avg), 





%%

function avg = get_avg(signal, idx, len)
avg = nan(length(idx), len);
T = length(signal);
for i = 1:length(idx)
    if idx(i) > T-len
        disp(idx(i))
        avg(i, 1:T-idx(i)+1) = signal(idx(i):T);
    else
        avg(i,:) = signal(idx(i):idx(i)+99);
    end
end
avg = mean(avg, 'omitnan');
end


function new_signal = addin(old_signal, func, idx)
T = length(old_signal);
for i = 1:length(idx)
    if idx(i) > T-100
        old_signal(idx(i):end) = old_signal(idx(i):end) + func(1:T-idx(i)+1);
    else
        old_signal(idx(i):idx(i)+99) = old_signal(idx(i):idx(i)+99) + func;
    end 
end
new_signal = old_signal;
end

