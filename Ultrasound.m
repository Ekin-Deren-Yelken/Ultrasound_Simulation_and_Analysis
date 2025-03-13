clear; clc; close all;
disp("running FINAL");
%% Global Variables
f_sample = 100e3; % 250 kHz Sampling Rate
f_o = 2e6; % 2 MHz Carrier frequencyue
pulse_duration = 10; % 10 second pulse duration
theta = deg2rad(60); % Ultrasound device Angle
c_sound = 1540; % [m/s] speed of sound in tissue

% Healthy Patient
healthy_PSV = 1; % [m/s] peak systolic velocity
healthy_EDV = 0.3;

file = 'heart_blood_flow_velocity_theta60.csv'; % 'heart_blood_flow_velocity_correct.csv';%'hbv.csv';%
data = readtable(file);
time = data{:,1};
doppler_shift = data{:,2};
blood_velocity_GT = data{:,3};
N = length(doppler_shift);

%% Generate the Received Doppler Ultrasound Signal
phi = 2 * pi * cumsum(doppler_shift) / f_sample; % Compute phase shift by integrating f_doppler
received_signal = cos(2 * pi * doppler_shift .* time + phi); % Modulated signal

%% Simulate Noise 
% Parameters
A_clutter = 1.5;  % Clutter amplitude
f_clutter = 0.5; % Clutter frequency
A_elec = 0.3;  % Power line noise amplitude
f_elec = 60;   % Power line frequency
sigma = 0.15;   % Thermal noise standard deviation
speckle_factor = 0.4; % Speckle noise intensity


% Clutter Tissue Noisetruncate = 
s_clutter = A_clutter * cos(2 * pi * f_clutter * time); % Low-frequency clutter
noisy_signal = received_signal + s_clutter;

% Electric Noise
line_noise = A_elec * sin(2 * pi * f_elec * time);
noisy_signal = noisy_signal + line_noise;

% Thermal Noise
thermal_noise = sigma * randn(size(time)); % White Gaussian noise
noisy_signal = noisy_signal + thermal_noise;

% Speckle Noise
speckle_noise = 1 + speckle_factor * randn(size(time)); % Multiplicative noise
noisy_signal = noisy_signal .* speckle_noise;

% Plot Signals (Time Domain)
start_trunk = 10e3;
truncate = 15e3; % for clarity

figure;
subplot(2,1,1);
plot(time(start_trunk:truncate), received_signal(start_trunk:truncate), 'k');
title('Original Doppler Signal (Truncated for Clarity)');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(2,1,2);
plot(time(start_trunk:truncate), noisy_signal(start_trunk:truncate), 'r');
title('Noisy Doppler Signal (Truncated for Clarity)');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

sgtitle('Original vs. Noisy Signal');

%% Compute Two-Sided FFT of Received Signal
Y = fft(noisy_signal); % Compute FFT
Y_mag = abs(Y) / N; % Normalize by N
freq_axis = (-f_sample/2 : f_sample/N : f_sample/2 - f_sample/N); 
Y_mag_shifted = fftshift(Y_mag); % Shift FFT to center it at 0 Hz

%% Peak Detection in Two-Sided FFT

function [filtered_freqs, filtered_pks] = detectProminentPeaks(Y_mag_shifted, freq_axis, prominence, neighborhood_threshold)
    % Peak Detection with Neighborhood Filtering
    % Inputs:
    %   Y_mag_shifted - Magnitude spectrum (two-sided FFT)
    %   freq_axis - Frequency vector corresponding to the FFT
    %   prominence - Minimum prominence for peak detection
    %   neighborhood_threshold - Frequency window for merging peaks (Hz)
    % Outputs:
    %   filtered_freqs - Prominent peak frequencies (Hz)
    %   filtered_pks - Corresponding magnitudes of the peaks

    % Find peaks in the magnitude spectrum
    [pks, locs] = findpeaks(Y_mag_shifted, 'MinPeakProminence', prominence); 
    prominent_freqs = freq_axis(locs); % Map peak indices to actual frequencies

    % Sort peaks by frequency
    [sorted_freqs, sort_idx] = sort(prominent_freqs);
    sorted_pks = pks(sort_idx);

    % Initialize filtered peak lists
    filtered_freqs = [];
    filtered_pks = [];

    i = 1;
    while i <= length(sorted_freqs)
        % Find all peaks within the ±neighborhood_threshold Hz range
        in_range_idx = find(abs(sorted_freqs - sorted_freqs(i)) <= neighborhood_threshold);

        % Get the strongest peak in this group
        [max_peak, max_idx] = max(sorted_pks(in_range_idx));
        filtered_freqs = [filtered_freqs; sorted_freqs(in_range_idx(max_idx))]; %#ok<AGROW>
        filtered_pks = [filtered_pks; max_peak]; %#ok<AGROW>

        % Skip to the next non-overlapping peak
        i = in_range_idx(end) + 1;
    end
end

[filtered_freqs_before, filtered_pks_before] = detectProminentPeaks(Y_mag_shifted, freq_axis, 0.1, 5);

% Plot Two-Sided FFT with Filtered Peaks
figure;
plot(freq_axis, 20*log10(Y_mag_shifted), 'k', 'LineWidth', 1.5); hold on;
scatter(filtered_freqs_before, 20*log10(filtered_pks_before), 'ro', 'DisplayName', 'Filtered Peaks');
title('Two-Sided FFT with Peak Detection');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend();
grid on;

%% High-Pass Filter to Remove Clutter
cutoff_freq = 50;  % Set cutoff at 50 Hz (can adjust to 100 Hz if needed)
filter_order = 4;  % Higher order = sharper cutoff
[b, a] = butter(filter_order, cutoff_freq / (f_sample/2), 'high');
filtered_signal = filtfilt(b, a, noisy_signal);

wall_filter_result = filtered_signal*1;

%% Notch Filter to Remove 60 Hz Power Line Noise 
notch_freq = 60;  % Target frequency (Hz)
notch_bandwidth = 1;  % Bandwidth of the notch filter (adjust if needed)

% Compute the lower and upper cutoff frequencies for the notch
low_cutoff = notch_freq - notch_bandwidth/2;
high_cutoff = notch_freq + notch_bandwidth/2;

% Normalize to Nyquist frequency
low_cutoff_norm = low_cutoff / (f_sample / 2);
high_cutoff_norm = high_cutoff / (f_sample / 2);

% Design a second-order Butterworth band-stop filter
[b_notch, a_notch] = butter(2, [low_cutoff_norm, high_cutoff_norm], 'stop');

% Apply Notch Filter to the High-Pass Filtered Signal
notch_filtered_signal = filtfilt(b_notch, a_notch, filtered_signal);
notch_filter_results = notch_filtered_signal*1;

%% Reconstruction via Continuous Wavelet Transform (CWT)
% Hilbert Transform Reconstruction with Constant Amplitude
% Assumes: notch_filtered_signal, f_sample, and time are already defined.

signal = notch_filtered_signal;
fs = f_sample;

% --- Wavelet Denoising ---
wname = 'db4';
level = floor(log2(length(signal)) - 1);
denoised_wavelet = wdenoise(signal, level, 'DenoisingMethod', "UniversalThreshold");
denoised_wavelet = denoised_wavelet - mean(denoised_wavelet);
% Optionally, if additional smoothing is desired, you can use sgolayfilt:
denoised_signal = sgolayfilt(denoised_wavelet, 3, 21);

% --- Hilbert Transform ---
analytic_signal = hilbert(denoised_signal);
inst_phase_hilbert = unwrap(angle(analytic_signal));
% Ensure the phase is a column vector (for element-wise operations)
inst_phase_hilbert = inst_phase_hilbert(:);

% --- Reconstruction with Constant Amplitude ---
% Instead of using the variable amplitude, we force it to 1.
reconstructed_signal = cos(inst_phase_hilbert);

% --- Plot the Results ---
figure;
subplot(3,1,1);
plot(time(start_trunk:truncate), denoised_signal(start_trunk:truncate));
title('Denoised Signal (Hilbert)');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(3,1,2);
plot(time, inst_phase_hilbert);
title('Instantaneous Phase (Hilbert)');
xlabel('Time (s)'); ylabel('Phase (radians)'); grid on;

subplot(3,1,3);
plot(time(start_trunk:truncate), reconstructed_signal(start_trunk:truncate));
title('Reconstructed Signal with Constant Amplitude (Hilbert)');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;


figure;
subplot(2,1,1);
plot(time(start_trunk:truncate), received_signal(start_trunk:truncate), 'k');
title('Original Doppler Signal (Truncated for Clarity)');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(2,1,2);
plot(time(start_trunk:truncate), reconstructed_signal(start_trunk:truncate));
title('Reconstructed Signal with Constant Amplitude (Hilbert)');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

%% Filters Step-by-Step
figure;
subplot(4,1,1);
plot(time(start_trunk:truncate), noisy_signal(start_trunk:truncate), 'r'); hold on;
title('Noisy Doppler Signal (Before Wall Filter)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4,1,2);
plot(time(start_trunk:truncate), wall_filter_result(start_trunk:truncate), 'b'); hold on;
title('Filtered Doppler Signal (Clutter Removed)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4,1,3);
plot(time(start_trunk:truncate), notch_filter_results(start_trunk:truncate), 'b'); hold on;
title('Filtered Doppler Signal (Electric Noise Removed)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4,1,4);
plot(time(start_trunk:truncate), denoised_signal(start_trunk:truncate), 'b'); hold on;
title('Filtered Doppler Signal (Thermal Noise Removed)');
xlabel('Time (s)');
ylabel('Amplitude');

%{
%% Signal Reconstruction Using Peak Frequencies
% Compute Two-Sided FFT AFTER Wavelet Denoising
Y_wavelet_filtered = fft(reconstructed_signal);
Y_mag_wavelet_filtered = abs(Y_wavelet_filtered) / N; % Normalize
Y_mag_wavelet_filtered_shifted = fftshift(Y_mag_wavelet_filtered);

% Detect Peaks After Wavelet Denoising
[filtered_freqs_wavelet, filtered_pks_wavelet] = detectProminentPeaks(Y_mag_wavelet_filtered_shifted, freq_axis, 0.1, 5);

% Compute FFT phase spectrum for detected frequencies
Y_phase_wavelet_filtered = angle(Y_wavelet_filtered); % Extract phase spectrum
Y_phase_wavelet_filtered_shifted = fftshift(Y_phase_wavelet_filtered); % Shift to match frequency axis

% Find phase corresponding to detected frequencies
filtered_phases_wavelet = zeros(size(filtered_freqs_wavelet)); % Initialize

for i = 1:length(filtered_freqs_wavelet)
    % Find the index closest to the detected frequency
    [~, idx] = min(abs(freq_axis - filtered_freqs_wavelet(i)));
    filtered_phases_wavelet(i) = Y_phase_wavelet_filtered_shifted(idx); % Store phase
end

% Step 1: Match Reconstruction Duration
reconstruction_duration = length(time); % Ensure it matches original time vector

% Ensure Peaks are Computed
if isempty(filtered_freqs_wavelet) || isempty(filtered_pks_wavelet)
    error('filtered_freqs_wavelet is empty! Peak detection might have failed.');
end

% Step 2: Use Peaks Detected After Wavelet Denoising
reconstruction_freqs = filtered_freqs_wavelet;  % Selected peak frequencies
reconstruction_amps = filtered_pks_wavelet;  % Corresponding amplitudes
reconstruction_phases = phi_t;  % Use extracted instantaneous phase φ(t)

% Step 3: Generate Reconstructed Signal
reconstructed_signal = zeros(size(time)); % Initialize reconstruction

for i = 1:length(reconstruction_freqs)
    reconstructed_signal = reconstructed_signal + ...
        reconstruction_amps(i) * cos(2 * pi * reconstruction_freqs(i) * time + reconstruction_phases(i));
end

% Step 4: Normalize Reconstructed Signal
reconstructed_signal = reconstructed_signal / max(abs(reconstructed_signal)) * max(abs(denoised_signal));

% Step 5: Store Final Result
final_reconstruction_result = reconstructed_signal;

% Step 6: Plot Comparison (Original vs. Reconstructed)
figure;
subplot(3,1,1);
plot(time(start_trunk:truncate), denoised_signal(start_trunk:truncate));
title('Denoised Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(time(start_trunk:truncate), final_reconstruction_result(start_trunk:truncate));
title('Reconstructed Signal (Using Extracted Frequencies & Phase)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,3);
plot(time(start_trunk:truncate), received_signal(start_trunk:truncate), 'k');
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
%}
%% Compute and Plot Spectrograms

window_length = 256; % Window size for STFT
overlap = 128; % Overlap between windows
nfft = 512; % FFT points

figure;
subplot(1,3,1);
spectrogram(noisy_signal, window_length, overlap, nfft, f_sample, 'yaxis');
title('Spectrogram of Noisy Doppler Signal');

subplot(1,3,2);
spectrogram(denoised_signal, window_length, overlap, nfft, f_sample, 'yaxis');
title('Spectrogram of Denoised Doppler Signal');

subplot(1,3,3);
spectrogram(reconstructed_signal, window_length, overlap, nfft, f_sample, 'yaxis');
title('Spectrogram of Reconstructed Doppler Signal');
%{%}

%% Compare Reconstructed vs. Original Signal

% Correlation Coefficient
corr_coeff_signal = corrcoef(received_signal, reconstructed_signal);
disp(['Correlation Coefficient (Reconstructed vs. Theory): ', num2str(corr_coeff_signal(1,2))]);

corr_coeff_signal_2 = corrcoef(received_signal, noisy_signal);
disp(['Correlation Coefficient (Noise vs. Theory): ', num2str(corr_coeff_signal_2(1,2))]);

analytic_signal = hilbert(reconstructed_signal);
inst_phase_hilbert = unwrap(angle(analytic_signal));

inst_freq = diff(inst_phase_hilbert)*f_sample/(2*pi);
inst_freq = [inst_freq; inst_freq(end)]; 

blood_velocity = (inst_freq * c_sound) ./ (2 * f_o * cos(theta));

% Optionally, take the absolute value if you want just the speed:
blood_velocity = abs(blood_velocity);

% Plot the blood velocity over time
figure;
plot(time, blood_velocity./100, 'LineWidth',1.5); hold on;
title('Blood Velocity vs. Time');
xlabel('Time (s)');
ylabel('Velocity (cm/s)');
grid on;

% Compute Doppler Shift Directly
doppler_shift_extracted = diff(phi) / (2 * pi * (1 / f_sample));
time_velocity = time(1:end-1); % Adjust time vector

% Step 3: Convert Doppler Shift to Blood Velocity
blood_velocity_extracted = (doppler_shift_extracted * c_sound) / (2 * f_o * cos(theta));

% Step 4: Plot Blood Velocity Over Time
plot(time, blood_velocity_GT-0.3, 'g', 'LineWidth', 1.5);
plot(time_velocity, blood_velocity_extracted-0.3, 'r', 'LineWidth', 1.5);
title('Extracted Blood Velocity');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
grid on;

%% Validation
theory_signal = received_signal; % Your theoretical reference signal

% Ensure signals have matching length
min_length = min(length(theory_signal), length(reconstructed_signal));
theory_resized = theory_signal(1:min_length);
reconstructed_resized = reconstructed_signal(1:min_length);
noisy_resized = noisy_signal(1:min_length); % Ensure noisy signal is also matched

% Mean Absolute Error (MAE)
MAE_recon = mean(abs(reconstructed_resized - theory_resized));
MAE_noise = mean(abs(noisy_resized - theory_resized));

% Root Mean Square Error (RMSE)
RMSE_recon = sqrt(mean((reconstructed_resized - theory_resized).^2));
RMSE_noise = sqrt(mean((noisy_resized - theory_resized).^2));

% Signal-to-Noise Ratio (SNR)
signal_power = mean(theory_resized.^2);
noise_power_recon = mean((theory_resized - reconstructed_resized).^2);
noise_power_noise = mean((theory_resized - noisy_resized).^2);

SNR_recon = 10 * log10(signal_power / noise_power_recon);
SNR_noise = 10 * log10(signal_power / noise_power_noise);

%% Display Results
fprintf('Error Metrics:\n');
fprintf(' - Mean Absolute Error (MAE):\n   Noise vs. Theory: %.6f\n   Reconstructed vs. Theory: %.6f\n', MAE_noise, MAE_recon);
fprintf(' - Root Mean Square Error (RMSE):\n   Noise vs. Theory: %.6f\n   Reconstructed vs. Theory: %.6f\n', RMSE_noise, RMSE_recon);
fprintf(' - Signal-to-Noise Ratio (SNR) (dB):\n   Noise: %.2f dB\n   Reconstructed: %.2f dB\n', SNR_noise, SNR_recon);

% CWT vs Hilbert

fb = cwtfilterbank('Wavelet', 'amor', 'SamplingFrequency', fs, 'SignalLength', length(denoised_signal));
[cfs, freq] = cwt(denoised_signal, 'FilterBank', fb);

% Extract phase information from the CWT coefficients
cfs_phase = angle(cfs);

% For each time step, select the phase corresponding to the maximum magnitude coefficient
[~, max_idx] = max(abs(cfs), [], 1);
inst_phase_cwt = zeros(1, length(denoised_signal));
for k = 1:length(denoised_signal)
    inst_phase_cwt(k) = cfs_phase(max_idx(k), k);
end
inst_phase_cwt = unwrap(inst_phase_cwt)';  % Convert to column vector

%% Comparison: Plotting the Two Phase Estimates
figure;
plot(time, inst_phase_hilbert, 'b-', 'DisplayName', 'Hilbert Phase');
hold on;
plot(time, inst_phase_cwt, 'r--', 'DisplayName', 'CWT Phase');
title('Comparison of Instantaneous Phase: Hilbert vs. CWT');
xlabel('Time (s)');
ylabel('Phase (radians)');
legend;
grid on;

%% Plot the Phase Difference
phase_difference = inst_phase_hilbert - inst_phase_cwt;
figure;
plot(time, phase_difference, 'k-', 'LineWidth', 1.5);
title('Phase Difference (Hilbert - CWT)');
xlabel('Time (s)');
ylabel('Phase Difference (radians)');
grid on;


%{
% Find the dominant frequency index for each time step (max magnitude at each time)
[~, max_idx] = max(abs(cfs), [], 1); % Find max magnitude for each time step
phi_t = zeros(1, length(time)); % Initialize phase array

for i = 1:length(time)
    phi_t(i) = cfs_phase(max_idx(i), i); % Extract phase at dominant frequency
end

% Step 3: Plot the Phase \( \phi(t) \)
figure;
plot(time, unwrap(phi_t)); % Unwrap phase for continuity
xlabel('Time (s)');
ylabel('Phase \( \phi(t) \) (radians)');
title('Instantaneous Phase Over Time');
grid on;


%}