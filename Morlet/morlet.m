load ../sampleEEGdata

time = -1:1/EEG.srate:1;

num_wavelets      =  5;  % number of frequency bands
lowest_frequency  =   2;  % in Hz
highest_frequency = 30;  % in Hz

frequencies=linspace(lowest_frequency,highest_frequency,num_wavelets);
% note: the "linspace" function creates linearly spaced numbers between the first and second 
% inputs, with the number of steps corresponding to the third input. 
figure, plot(frequencies,'-*')
xlabel('Frequency order')
ylabel('Frequency in Hz')


% initialize wavelet family
wavelet_family = zeros(num_wavelets,length(time));
 
% Loop through frequencies and make a family of wavelets.
for fi=1:num_wavelets
    
    % create a sine wave at this frequency
    sinewave = exp(2*1i*pi*frequencies(fi).*time); % the "1i" makes it a complex wavelet
    
    % create a Gaussian window
    gaus_win = exp(-time.^2./(2*(6/(2*pi*frequencies(fi)))^2));
    
    % create wavelet via element-by-element multiplication of the sinewave and gaussian window
    wavelet_family(fi,:) = sinewave.*gaus_win;
end

% Plot a few wavelets
figure
subplot(2,1,1)
plot(time,real(wavelet_family(1:round(rand*30):end,:))') 
title('A few wavelets...')
 
% Note that in the subplot command you don't need commas if you have fewer than 10 
% rows/cols and if you are not using any variables.
subplot(212)
plot(time,real(wavelet_family(5,:)))
hold on
plot(time,imag(wavelet_family(5,:)),':')
title('Real and imaginary parts of one wavelet')
legend({'real';'imaginary'})


% finally, image the wavelet family.
figure
imagesc(time,frequencies,real(wavelet_family))
axis xy % equivalent to "set(gca,'ydir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')


%% Question 2

% EEG data from one trial (electrode FCz)
eegdata = squeeze(EEG.data(47,:,10));

wavelet_conv_data = zeroes(num_wavelets, );

% convolve with data
% compute Gaussian
for fi=1:num_wavelets
    n_conv = length(wavelet_family(fi, :)) + EEG.pnts - 1;
    fft_w = fft(wavelet_family(fi, :), n_conv);
    fft_e = fft(eegdata, n_conv);

    ift   = ifft(fft_e.*fft_w,n_conv)*sqrt(s)/10; % sqrt... is an empirical scaling factor that works here
    wavelet_conv_data(fi:) = real(ift(halfwaveletsize:end-halfwaveletsize+1));
end


