clear
close all
clc


m = 5;
fs = 1000;
t = -1.5:1/fs:1.5;

% frequencies which we are going to use
wavelet_family = linspace(2, 30, 5);

load sampleEEGdata.mat

data = EEG.data(:,:, 47);

% Preallocate array for wavelets
wavelets = zeros(length(wavelet_family), length(t));

labels = []; 

for wi = 1:length(wavelet_family)
    w = morlet_wavelet(wavelet_family(wi), m, t-mean(t));
    wavelets(wi, :) = w;
    figure(1);
    plot(t, real(w));
    hold on;
end

hold off;


function w = morlet_wavelet(f0, m, t)
    sigma = m / (2*pi*f0);
    w = exp(2*1i*pi*f0*t) .* exp(-t.^2/(2*sigma^2));
    w = w ./ sqrt(sum(abs(w).^2)) / sqrt(sigma); % normalize energy and scale
end