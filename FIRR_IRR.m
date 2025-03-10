% FIR & IIR Filtering Example
% Designs and applies FIR & IIR filters on a noisy signal

clc; clear; close all;

% Generate a noisy sine wave
fs = 1000; % Sampling frequency (Hz)
t = 0:1/fs:1; % 1 second time vector
signal_clean = sin(2*pi*10*t) + 0.5*sin(2*pi*50*t); % Clean signal
noise = 0.3 * randn(size(signal_clean)); 
signal_noisy = signal_clean + noise; % Noisy signal

% FIR Filter Design (Window Method)
fir_order = 50;
cutoff_fir = 30/(fs/2); % Normalize cutoff frequency
fir_coeffs = fir1(fir_order, cutoff_fir, 'low'); % Low-pass FIR filter
signal_fir = filtfilt(fir_coeffs, 1, signal_noisy);

% IIR Filter Design (Butterworth)
iir_order = 4;
cutoff_iir = 30/(fs/2); % Normalize cutoff frequency
[b, a] = butter(iir_order, cutoff_iir, 'low'); % Low-pass IIR filter
signal_iir = filtfilt(b, a, signal_noisy);

% Plot results
figure;
subplot(3,1,1);
plot(t, signal_noisy, 'r');
title('Noisy Signal');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(3,1,2);
plot(t, signal_fir, 'b');
title('FIR Filtered Signal');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(3,1,3);
plot(t, signal_iir, 'g');
title('IIR Filtered Signal');
xlabel('Time (s)'); ylabel('Amplitude');

sgtitle('Comparison of FIR and IIR Filtering');

% Comments on filter performance:
% - FIR filters have a linear phase response but require a high order.
% - IIR filters are more efficient (lower order), but can introduce phase distortion.
