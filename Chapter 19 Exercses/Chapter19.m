clear;

load sampleEEGdata.mat

channels = ["P1", "P3", "P7"];

[num_channels, num_timepoints, num_trials] = size(EEG.data);
power = zeros(3, num_trials, num_timepoints);

indices = [];

for i=3
    indices = [indices, EEG.chanlocs.urchan(find(EEG.chanlocs.labels==channels(:, i)))];
end