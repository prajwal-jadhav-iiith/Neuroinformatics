% Question 2

% step size
dt = 1 / 1000;

% time vector with step size
time = (0:dt:4);

% below we create 4 waves with frequencies 2, 8, 12, 25 Hz
freq1 = 2;
sine_wave1 = 1 * sin(2 * pi * freq1 * time + pi / 4);

freq2 = 8;
sine_wave2 = 2 * sin(2 * pi * freq2 * time + pi / 3);

freq3 = 12;
sine_wave3 = 3 * sin(2 * pi * freq3 * time + pi / 2);

freq4 = 25;
sine_wave4 = 4 * sin(2 * pi * freq4 * time + pi);

% average all the waves to create a new wave
mixed_wave = (sine_wave1 + sine_wave2 + sine_wave3 + sine_wave4) / 4;

% plot a 5*1 figure
figure(1);
subplot(5, 1, 1);
plot(time, sine_wave1);
title('Sine Wave 1: 2 Hz');
xlabel('Time');
ylabel('Amplitude');

subplot(5, 1, 2);
plot(time, sine_wave2);
title('Sine Wave 2: 8 Hz');
xlabel('Time');
ylabel('Amplitude');

subplot(5, 1, 3);
plot(time, sine_wave3);
title('Sine Wave 3: 12 Hz');
xlabel('Time');
ylabel('Amplitude');

subplot(5, 1, 4);
plot(time, sine_wave4);
title('Sine Wave 4: 25 Hz');
xlabel('Time');
ylabel('Amplitude');

subplot(5, 1, 5);
plot(time, mixed_wave);
title('Mixed Wave');
xlabel('Time');
ylabel('Amplitude');

% run the mixed wave through fast fourier transform
fft_transform = fft(mixed_wave);

% calculate the length of the signal, used to normalize
length_of_signal = length(mixed_wave);

% normalize the fft output with the length of the signal and
% extract the real part of complex numbers by using abs() function
first_half = abs(fft_transform / length_of_signal);

% since the first_half variable contains symmetric data we can
% select only first half of the data to focus on positive frequencies
first_half = first_half(1:length_of_signal/2 + 1);

% multiply by 2 to compensate for the loss of data in the above step
% ignore first and last elements
first_half(2:end-1) = 2*first_half(2:end-1);

% Create a frequency vector for a signal based on its length
frequency_vector = 1000*(0:(length_of_signal / 2)) / length_of_signal;

% plot the figure
figure(2);
plot(frequency_vector,first_half);
title('Mixed wave decoded in its constituents')
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% Question 3: Add random white noise and then recover

% choose a random noise amplitude
noise_amplitude = 0.5;

% Create the mixed wave with added noise
noise = noise_amplitude * randn(size(mixed_wave));
mixed_wave_with_noise = mixed_wave + noise;

% plot the figure
figure(3);
plot(time, mixed_wave_with_noise);
title('Mixed wave with noise');
xlabel('Time');
ylabel('Amplitude');

fft_transform_with_noise = fft(mixed_wave_with_noise);

first_half_with_noise = abs(fft_transform_with_noise / length_of_signal);
first_half_with_noise = first_half_with_noise(1:length_of_signal/2 + 1);
first_half_with_noise(2:end-1) = 2*first_half_with_noise(2:end-1);

figure(4);
plot(frequency_vector, first_half_with_noise);
title('Mixed wave with noise decoded in its constituents');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

