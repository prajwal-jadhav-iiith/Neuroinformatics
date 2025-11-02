load sampleEEGdata.mat

EEG = EEG;

chan2use = 'fcz';
min_freq = 3;
max_freq = 30;
num_frex = 20;

chanindex = find(strcmpi({EEG.chanlocs.labels}, chan2use));

frex = logspace(log10(min_freq), log10(max_freq), num_frex);

for ev = 1:size(EEG.data, 3)
    [cfs(:, :, ev), ~] = morletWaveletTransform(EEG.data(chanindex, : , ev), EEG.srate, frex, 6, 2);
end

power_all = abs(cfs).^2;

% Edge Trimming
time_s = dsearchn(EEG.times', -500);
time_e = dsearchn(EEG.times', 1200);

eegpower = power_all(:, time_s:time_e, :);
tftimes = EEG.times(time_s:time_e);
nTimepoints = length(tftimes);

baseidx = [dsearchn(tftimes', -500), dsearchn(tftimes', -100)];
% do trial level baseline normalization, store it in realmean
realbaseline = squeeze(mean(eegpower(:, baseidx(1):baseidx(2), :), 2));

realmean = zeros(20,436,99);
%% 

% apply baseline normalization on eegpower
for i = (1:20)
    for j = (1:99)
        realmean(i, :, j) = eegpower(i, :, j) - realbaseline(i, j);
    end
end

trial_average = mean(realmean, 3);

permuted_trials = zeros(size(trial_average, 1), size(trial_average, 2), 1000);
for k = 1:1000
    permuted_trials(:, :, k) = circshift(trial_average, randi(size(trial_average, 2)), 2);
end

mean_permuted = mean(permuted_trials, 3);
std_permuted = std(permuted_trials, 0, 3);

z_score = (realmean - mean_permuted) ./ std_permuted;
