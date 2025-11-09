load sampleEEGdata.mat

disp(EEG);

sampling_freq = EEG.srate;

theta_vector = [4 8] / (sampling_freq / 2);

[b, a] = butter(4, theta_vector, "bandpass");

data = EEG.data;

permuted_data = permute(double(data), [2, 1, 3]);

theta_permuted = filtfilt(b, a, permuted_data);

theta_filtered = permute(theta_permuted, [2, 1, 3]);

theta_power = abs(hilbert(theta_filtered)).^2;

base_idx = dsearchn(EEG.times', [-500 -100]');

baseline_power = mean(mean(theta_power(:, base_idx(1):base_idx(2), :), 2), 3);

dB = 10 * log10(theta_power ./ baseline_power);