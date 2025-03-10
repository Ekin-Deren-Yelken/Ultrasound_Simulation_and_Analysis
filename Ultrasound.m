% Clear Workspace and Load Data
clear; clc; close all;

%% Parameters
f_sample = 1e6;  % Sampling frequency (1 MHz)
f_transmit = 500e3;  % Frequency of transmitted ultrasound wave (500 kHz)
v = 0.5;    % Blood flow velocity (m/s)
c = 1540;   % Speed of sound in tissue (m/s)
theta = deg2rad(45);  % Angle in radians
duration = 1e-2; % Signal duration (10 ms)
time = linspace(0, duration, f_sample * duration);

%% Generate Original Doppler Signal
v_pulsatile = v + 0.2 * sin(2 * pi * (75/60) * time);  % Simulated heartbeat pulsations
f_doppler = (2 * f_transmit * v_pulsatile .* cos(theta)) / c;  % Doppler shift
f_received = f_transmit + f_doppler;
received_signal = sin(2 * pi * (f_received .* time));  % Original Doppler signal (no noise)

%% Add Noise (Realistic Simulation)
noise_level = 0.3; % Gaussian noise
gaussian_noise = noise_level * randn(size(received_signal));

% Powerline Noise (50 Hz)
line_freq = 50;  
A = 0.4;  
line_noise = A * sin(2 * pi * line_freq * time) +  + 0.5*sin(2*pi*50*time);

% Motion Artifact Noise (0.5 Hz slow variations)
motion_artifact = A * sin(2 * pi * 0.5 * time);

% Harmonic Distortion
harmonic_amplitude = 0.1;
harmonic_noise = harmonic_amplitude * (sin(2 * pi * 2 * f_transmit * time) + ...
                                       sin(2 * pi * 3 * f_transmit * time));

% Add all noise components
received_signal_noisy = received_signal + gaussian_noise + line_noise + motion_artifact + harmonic_noise;

% Normalize noisy signal
received_signal_noisy = ( received_signal_noisy / max(abs(received_signal_noisy)) ); % - mean(received_signal_noisy);


%% Plot Original vs. Noisy Signal
figure;
subplot(2,1,1);
plot(time, received_signal, 'k');
title('Original Doppler Signal');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(2,1,2);
plot(time, received_signal_noisy, 'r');
title('Noisy Doppler Signal');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

sgtitle('Original vs. Noisy Signal');

%% Wavelet Denoising (Improved)
waveletName = 'db4';  
level = 7;  

% Perform wavelet decomposition
[coeffs, levels] = wavedec(received_signal_noisy, level, waveletName);

threshold = median(abs(coeffs)) / 0.6745;  
coeffs = wthresh(coeffs, 's', threshold);  % Apply soft thresholding
received_signal_denoised = waverec(coeffs, levels, waveletName);  % Reconstruct


%% FIR Bandpass Filter (Redesigned)
% Instead of a low-pass filter, use a **bandpass filter** to preserve Doppler shifts
filter_order = 60;  % FIR filter order
nyquist_freq = f_sample / 2;

% Ensure bandpass filter cutoff frequencies are within the valid range (0 < freq < 1)
bandpass_low = max(1e3, f_transmit - 20e3) / nyquist_freq;  % Ensure positive
bandpass_high = min(nyquist_freq - 1e3, f_transmit + 20e3) / nyquist_freq;  %

% Design FIR Bandpass Filter
fir_coeffs = fir1(filter_order, [bandpass_low bandpass_high], 'bandpass');

% Apply FIR filtering
filtered_signal = filter(fir_coeffs, 1, received_signal_denoised);

%% Frequency Analysis using FFT
L = length(received_signal_noisy);  % Use actual signal length
freq_axis = linspace(0, f_sample, L);  % Frequency axis from 0 to f_sample

fft_original = abs(fft(received_signal));  % FFT of original signal
fft_noisy = abs(fft(received_signal_noisy));  % FFT before denoising
fft_denoised = abs(fft(received_signal_denoised));  % FFT after wavelet denoising
fft_filtered = abs(fft(filtered_signal));  % FFT after FIR filtering

% Plot FFTs
figure;

subplot(4,1,1);
plot(freq_axis, abs(fft_original), 'k');
title('FFT of Original Signal');
xlabel('Frequency (Hz)'); ylabel('Magnitude'); grid on;

subplot(4,1,2);
plot(freq_axis, abs(fft_noisy), 'r');
title('FFT of Noisy Signal');
xlabel('Frequency (Hz)'); ylabel('Magnitude'); grid on;

subplot(4,1,3);
plot(freq_axis, abs(fft_denoised), 'b');
title('FFT After Wavelet Denoising');
xlabel('Frequency (Hz)'); ylabel('Magnitude'); grid on;

subplot(4,1,4);
plot(freq_axis, abs(fft_filtered), 'g');
title('FFT After FIR Bandpass Filtering');
xlabel('Frequency (Hz)'); ylabel('Magnitude'); grid on;

sgtitle('FFT Analysis: Original vs. Noisy vs. Denoised vs. Filtered');

%% Spectrogram to visualize frequency shifts
figure;
spectrogram(filtered_signal, 256, 200, 512, f_sample, 'yaxis');
title('Spectrogram of Filtered Doppler Signal');
xlabel('Time (s)'); ylabel('Frequency (Hz)');
colorbar;

%% SNR Comparison
snr_before = snr(received_signal, received_signal_noisy - received_signal);
snr_after = snr(received_signal, filtered_signal - received_signal);
fprintf('SNR Before Denoising: %.2f dB\n', snr_before);
fprintf('SNR After Denoising: %.2f dB\n', snr_after);

%% Plot Comparison: Noisy vs. Wavelet Denoised vs. FIR Filtered
figure;
subplot(3,1,1);
plot(time, received_signal_noisy, 'r');
title('Noisy Signal');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(3,1,2);
plot(time, received_signal_denoised, 'b');
title('Wavelet Denoised Signal');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(3,1,3);
plot(time, filtered_signal, 'g');
title('Final FIR Filtered Signal');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

sgtitle('Comparison of Noisy vs. Denoised vs. FIR Filtered Signal');

%% Correlation Analysis
corr_noisy = corr(received_signal', received_signal_noisy');  % Original vs Noisy
corr_denoised = corr(received_signal', received_signal_denoised');  % Original vs Wavelet Denoised
corr_filtered = corr(received_signal', filtered_signal');  % Original vs FIR Filtered

% Display Correlation Results
fprintf('Correlation (Original vs Noisy): %.4f\n', corr_noisy);
fprintf('Correlation (Original vs Wavelet Denoised): %.4f\n', corr_denoised);
fprintf('Correlation (Original vs FIR Filtered): %.4f\n', corr_filtered);

%% Introduce Phase Distortion
% Convert noisy signal to frequency domain
fft_noisy = fft(received_signal_noisy);

% Create a phase shift (example: 90-degree shift at all frequencies)
phase_shift = exp(1j * pi/2 * (0:length(fft_noisy)-1));

% Apply the phase shift
fft_noisy_shifted = fft_noisy .* phase_shift;

% Convert back to time domain
received_signal_phase_distorted = real(ifft(fft_noisy_shifted));

%% Process the Phase-Distorted Signal with Wavelet & FIR Filter
% Apply Wavelet Denoising
received_signal_denoised = wdenoise(received_signal_phase_distorted, 3, ...
    'Wavelet', 'db4', 'DenoisingMethod', 'UniversalThreshold');

% Design a Linear Phase FIR Bandpass Filter
filter_order = 40;
bandpass_low = max(1e3, f_transmit - 20e3) / nyquist_freq;  % Ensure positive
bandpass_high = min(nyquist_freq - 1e3, f_transmit + 20e3) / nyquist_freq;  %
fir_coeffs = fir1(filter_order, [bandpass_low bandpass_high], 'bandpass');

% Apply FIR filter
filtered_signal = filtfilt(fir_coeffs, 1, received_signal_denoised);

%% Plot Signals to Compare Phase Distortion & Correction
figure;

subplot(4,1,1);
plot(time, received_signal_noisy, 'r');
title('Noisy Signal (Before Phase Distortion)');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(4,1,2);
plot(time, received_signal_phase_distorted, 'm');
title('Phase-Distorted Signal');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(4,1,3);
plot(time, received_signal_denoised, 'b');
title('Wavelet Denoised Signal (No Phase Correction)');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(4,1,4);
plot(time, filtered_signal, 'g');
title('Final FIR Filtered Signal (Phase Correction Applied)');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

sgtitle('Comparison of Phase Distortion & Correction');

%% Check Phase Response of FIR Filter
figure;
freqz(fir_coeffs, 1, 1024, f_sample);
title('Frequency & Phase Response of FIR Filter');
