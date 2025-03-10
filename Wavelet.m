% Wavelet Denoising Example (Fixed Version)
% Generates a synthetic signal with noise and applies wavelet denoising

clc; clear; close all;

% Generate a clean signal (sum of two sine waves)
fs = 1000; % Sampling frequency (Hz)
t = 0:1/fs:1; % Time vector (1 second)
signal_clean = sin(2*pi*10*t) + 0.5*sin(2*pi*50*t); 

% Add Gaussian noise
noise = 0.3 * randn(size(signal_clean));
signal_noisy = signal_clean + noise;

% Perform wavelet denoising using Daubechies wavelets
waveletName = 'db4'; % Daubechies wavelet


% Use wdenoise for automatic threshold selection and denoising
level_mid = 2; % Decomposition level
signal_denoised = wdenoise(signal_noisy, level_mid, ...
    'Wavelet', waveletName, 'DenoisingMethod', 'UniversalThreshold');

% Apply low-level wavelet denoising (Level 1)
level_low = 1;
signal_denoised_low = wdenoise(signal_noisy, level_low, ...
    'Wavelet', waveletName, 'DenoisingMethod', 'UniversalThreshold');

% Apply high-level wavelet denoising (Level 4)
level_high = 4;
signal_denoised_high = wdenoise(signal_noisy, level_high, ...
    'Wavelet', waveletName, 'DenoisingMethod', 'UniversalThreshold');

% Plot results
figure;

subplot(5,1,1);
plot(t, signal_clean, 'k', 'LineWidth', 1.2);
title('Original Signal');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(5,1,2);
plot(t, signal_noisy, 'r');
title('Noisy Signal');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(5,1,3);
plot(t, signal_denoised_low, 'b');
title(['Denoised Signal - Low Decomposition Level (', num2str(level_low), ')']);
xlabel('Time (s)'); ylabel('Amplitude');

subplot(5,1,4);
plot(t, signal_denoised_high, 'g');
title(['Denoised Signal - High Decomposition Level (', num2str(level_high), ')']);
xlabel('Time (s)'); ylabel('Amplitude');

subplot(5,1,5);
plot(t, signal_denoised, 'black');
title(['Denoised Signal - Mid Decomposition Level (', num2str(level_mid), ')']);
xlabel('Time (s)'); ylabel('Amplitude');

sgtitle('Wavelet Denoising: Low vs. High Decomposition Levels');