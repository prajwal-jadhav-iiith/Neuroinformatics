clear
close all
clc


m = 5;
fs = 1000;
t = 0:1/fs:100;

wavelet_family = 2:5:30;

load ../sampleEEGdata.mat

data = EEG.data(:,:, 47);

% Preallocate array for wavelets
wavelets = zeros(length(wavelet_family), length(t));

for wi = 1:length(wavelet_family)
    % wf = wavelet_family(wi);
    % sigma = m/(2*pi*wf);
    % tw = -6*sigma : 1/fs : 6*sigma;
    % w = morlet_wavelet(wf, m, tw);
    % 
    % figure(wi)
    % plot(t, real(w))

    w = morlet_wavelet(wavelet_family(wi), m, t-mean(t));

    figure(wi);

    plot(t, real(w));

    
end


function w = morlet_wavelet(f0, m, t)
    sigma = m / (2*pi*f0);
    w = exp(2*1i*pi*f0*t) .* exp(-t.^2/(2*sigma^2));
    w = w ./ sqrt(sum(abs(w).^2)) / sqrt(sigma); % normalize energy and scale
end