% Question 1
kernel_size = 15;
domain = linspace(-1, 1, kernel_size);

inverted_u_kernel = 1 - domain.^2;

inverted_u_kernel = inverted_u_kernel / sum(inverted_u_kernel);

disp(inverted_u_kernel);

figure(1);
plot(inverted_u_kernel, '-o', 'LineWidth', 2);
title('Inverted U Kernel');
grid on;

decay_rate = 0.5;

x = 0:(kernel_size-1);

decay_kernel = exp(-decay_rate * x);

decay_kernel = decay_kernel / sum(decay_kernel);

disp(decay_kernel);

figure(2);
plot(decay_kernel, '-o', 'LineWidth', 2);
title('1D Exponential Decay Kernel');
xlabel('Kernel Index');
ylabel('Weight');
grid on;

eeg = load("sampleEEGdata.mat").EEG;
data = eeg.data(47, 1:50, 1);

data_len = length(data);
conv_len = data_len + kernel_size - 1;

convolution_result_u = zeros(1,conv_length);

for i = 1:conv_len
    for j = 1:kernel_size
        if (i - j + 1 > 0) && (i - j + 1 <= data_len)
            convolution_result_u(i) = convolution_result_u(i) + data(i - j + 1) * inverted_u_kernel(j);
        end
    end
end


convolution_result_decay = zeros(1, conv_len);

for i = 1:conv_len
    for j = 1:kernel_size
        if (i - j + 1 > 0) && (i - j + 1 <= data_len)
            convolution_result_decay(i) = convolution_result_decay(i) + data(i - j + 1) * decay_kernel(j);
        end
    end
end


% Plot results for the Inverted U Kernel
figure(3);
sgtitle('Convolution with Inverted U Kernel');

subplot(3, 1, 1);
plot(inverted_u_kernel, '-o', 'LineWidth', 1.5);
title('Inverted U Kernel');
xlabel('Kernel Index');
ylabel('Weight');
grid on;

subplot(3, 1, 2);
plot(data, '-o', 'LineWidth', 1.5);
title('EEG Data (Fcz, 50 time-points)');
xlabel('Time Point');
ylabel('Amplitude (\muV)');
grid on;

subplot(3, 1, 3);
plot(convolution_result_u, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
title('Convolution Result');
xlabel('Time Point');
ylabel('Convolved Signal');
grid on;


% Plot results for the Exponential Decay Kernel
figure(4);
sgtitle('Convolution with Exponential Decay Kernel');

subplot(3, 1, 1);
plot(decay_kernel, '-o', 'LineWidth', 1.5);
title('Exponential Decay Kernel');
xlabel('Kernel Index');
ylabel('Weight');
grid on;

subplot(3, 1, 2);
plot(data, '-o', 'LineWidth', 1.5);
title('EEG Data (Fcz, 50 time-points)');
xlabel('Time Point');
ylabel('Amplitude (\muV)');
grid on;

subplot(3, 1, 3);
plot(convolution_result_decay, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
title('Convolution Result');
xlabel('Time Point');
ylabel('Convolved Signal');
grid on;

data = eeg.data(47, :, 1);

fft_len = data_len + kernel_size - 1;

data_fft_u = fft(data, fft_len);
kernel_fft_u = fft(inverted_u_kernel, fft_len);

convolution_fft_u = data_fft_u .* kernel_fft_u;
convolution_result_u = ifft(convolution_fft_u, fft_len);

data_fft_decay = fft(data, fft_len);
kernel_fft_decay = fft(decay_kernel, fft_len);
convolution_fft_decay = data_fft_decay .* kernel_fft_decay;
convolution_result_decay = ifft(convolution_fft_decay, fft_len);

convolution_result_u = real(convolution_result_u);
convolution_result_decay = real(convolution_result_decay);


figure(5);
sgtitle('FFT Convolution with Inverted U Kernel');

subplot(2, 1, 1);
plot(eeg.times, data);
title('Original EEG Data (Full Duration)');
xlabel('Time (ms)');
ylabel('Amplitude (\muV)');
grid on;

subplot(2, 1, 2);
conv_times = linspace(eeg.times(1), eeg.times(end), fft_len);
plot(conv_times, convolution_result_u);
title('Convolution Result');
xlabel('Time (ms)');
ylabel('Convolved Signal');
grid on;
xlim(eeg.times([1 end]));

figure(6);
sgtitle('FFT Convolution with Decay');

subplot(2, 1, 1);
plot(eeg.times, data);
title('Original EEG Data (Full Duration)');
xlabel('Time (ms)');
ylabel('Amplitude (\muV)');
grid on;

subplot(2, 1, 2);
conv_times = linspace(eeg.times(1), eeg.times(end), fft_len);
plot(conv_times, convolution_result_decay);
title('Convolution Result');
xlabel('Time (ms)');
ylabel('Convolved Signal');
grid on;
xlim(eeg.times([1 end]));