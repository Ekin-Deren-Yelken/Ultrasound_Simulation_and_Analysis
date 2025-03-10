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
level = 4; % Decomposition level

% Use wdenoise for automatic threshold selection and denoising
signal_denoised = wdenoise(signal_noisy, level, ...
    'Wavelet', waveletName, 'DenoisingMethod', 'UniversalThreshold');

% Plot results
figure;
subplot(3,1,1);
plot(t, signal_clean, 'k', 'LineWidth', 1.2);
title('Original Signal');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(3,1,2);
plot(t, signal_noisy, 'r');
title('Noisy Signal');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(3,1,3);
plot(t, signal_denoised, 'b');
title('Denoised Signal (Wavelet)');
xlabel('Time (s)'); ylabel('Amplitude');

sgtitle('Wavelet Denoising of Noisy Signal');

