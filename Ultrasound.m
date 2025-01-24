% Clear Workspace and Load Data
clear; clc; close all;

% Parameters
f_sample = 10e6;  % Sampling frequency (10 MHz)
f_transmit = 2e6;  % Frequency of wave transmitted by the transducer (2 MHz)
v = 0.5;    % Blood flow velocity (m/s)
c = 1540;   % Speed of sound in tissue (m/s)
theta = deg2rad(45);  % Angle in radians

duration = 1e-3; % duration in seconds

fprintf(['Doppler Shift Analysis Parameters:\n' ...
    'Sampling frequency = %.2f [Hz]\n' ...
    'Transmitted frequency = %.2f [Hz]\n' ...
    'Blood flow velocity = %.2f [m/s]\n' ...
    'Speed of sound in tissue = %.2f [m/s]\n' ...
    'Transmission Angle = %.2f [radians]\n' ...
    'Duration = %.4f [sec]\n'], f_sample, f_transmit,v, c, theta, duration);

% Time vector
time = linspace(0, duration, f_sample * duration);  % 1 ms duration

%% Simulate pulsatile blood flow variation due to heartbeat
v_base = 0.5;       % Base blood flow velocity (m/s)
v_pulse_amplitude = 0.2;  % Velocity fluctuation due to heartbeat

heart_rate = 75 / 60;  % 75 BPM in Hz (1.25 Hz)
v_pulsatile = v_base + v_pulse_amplitude * sin(2 * pi * heart_rate * time);

% Doppler shift calculation
f_doppler = (2 * f_transmit * v_pulsatile .* cos(theta)) / c;
fprintf('\n\nDoppler Shift Frequency (Mean): %.2f Hz\n', mean(f_doppler));

f_recieved = f_transmit + f_doppler;

alpha = 0.33;

% Modelling Signal with Doppler Shift
received_signal = alpha*sin(2 * pi * (f_recieved .* time));
transmitted_signal = alpha*sin(2 * pi * (f_transmit .* time));

% Add noise to simulate real-world conditions

% Gaussian Noise
noise_level = 0.2; % Noise amplitude
gaussian_noise = noise_level * randn(size(received_signal));
% Powerline Noise
line_freq = 50;  % Powerline frequency (50 Hz or 60 Hz)
A = 0.5;
line_noise = A * sin(2 * pi * line_freq * time);  % Low amplitude noise
% Motion Artifact Noise
motion_freq = 0.5; % Hz
motion_artifact = A * sin(2 * pi * motion_freq * time);  % Slow movement artifact
% Simulate Harmonic Distortion
harmonic_amplitude = 0.05; 
harmonic_noise = harmonic_amplitude * (sin(2 * pi * 2 * f_transmit * time) + ...
                                       sin(2 * pi * 3 * f_transmit * time) + ...
                                       sin(2 * pi * 4 * f_transmit * time) + ...
                                       sin(2 * pi * 5 * f_transmit * time) + ...
                                       sin(2 * pi * 6 * f_transmit * time));
% Simulate frequency drift over time
drift_rate = 500;  % Frequency drift rate in Hz
f_drift = f_transmit + drift_rate * time / max(time);
drift_signal = sin(2 * pi * f_drift .* time);

% Recieved Signal
received_signal_noisy = received_signal + gaussian_noise + line_noise + motion_artifact + harmonic_noise +  0.1 * drift_signal;

% Normalize signal
received_signal_noisy = received_signal_noisy / max(abs(received_signal_noisy));

L = 6000;

% Plot the signal
figure;
plot(time(1:L), transmitted_signal(1:L));
title('Transmitted Signal');
xlabel('Time (s)');
ylabel('Amplitude');

figure;
plot(time(1:L), received_signal_noisy(1:L), 'r');
title('Recieved Signal');
xlabel('Time (s)');
ylabel('Amplitude');

%% Frequency analysis using FFT
N = length(received_signal_noisy);
window = hamming(length(received_signal_noisy))';
frequency_content = fft(received_signal_noisy );%.* window);
freqs = (0:N-1)*(f_sample/N);
figure;
plot(freqs(1:end/2), abs(frequency_content(1:end/2)));
title('FFT Analysis of Doppler Shift');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Mark the detected Doppler frequency
hold on;
[max_amp, idx] = max(abs(frequency_content(1:N/2)));  % Find peak frequency
detected_freq = freqs(idx);  % Get frequency value corresponding to peak
stem(detected_freq, max_amp, 'r');
legend('FFT Spectrum', 'Doppler Shift');

%% Spectrogram to visualize frequency shift over time
figure;
spectrogram(received_signal_noisy, 256, 200, 512, f_sample, 'yaxis');
title('Spectrogram of Doppler Shifted Signal');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

%% Filtering Techniques to Isolate Blood Flow

% Notch filter parameters
notch_freq = 50;   % Powerline frequency (Hz)
Q_factor = 35;     % Quality factor (adjust to control sharpness)
w0 = (2 * pi * notch_freq) / f_sample;
bw = w0 / 35;              % Bandwidth
if bw < 1e-4
    bw = 1e-3;  % Set a minimum bandwidth to prevent instability
end

% Compute filter coefficients for notch filter
alpha = sin(w0) / (2 * Q_factor);
b = [1, -2*cos(w0), 1];  % Numerator coefficients
a = [1 + alpha, -2*cos(w0), 1 - alpha];  % Denominator coefficients

% Apply the filter
filtered_signal_notch = filtfilt(b, a, received_signal_noisy);

% Design a high-pass filter to remove motion artifacts (<5 Hz)
cutoff_freq_hp = 5;  % Cutoff frequency at 5 Hz
[b_hp, a_hp] = butter(2, cutoff_freq_hp / (f_sample / 2), 'high');

% Apply high-pass filter to remove low-frequency motion artifacts
filtered_signal_hp = filter(b_hp, a_hp, filtered_signal_notch);

% Design a low-pass filter to remove high-frequency noise (>5000 Hz)
cutoff_freq_lp = f_transmit + 3000;  % Cutoff frequency slightly above Doppler frequency
[b_lp, a_lp] = butter(4, cutoff_freq_lp / (f_sample / 2), 'low');

% Apply low-pass filter to remove high-frequency noise
filtered_signal_final = filter(b_lp, a_lp, filtered_signal_hp);

% Plot comparison
figure;
plot(time(1:L), filtered_signal_final(1:L));
title('After Low-Pass Filtering (Final Cleaned Signal)');
xlabel('Time (s)');
ylabel('Amplitude');

%%
% Perform FFT analysis of the filtered signal
N = length(filtered_signal_final);

% Plot frequency spectrum
frequency_content_filtered = fft(filtered_signal_final);% .* window);
freqs = (0:N-1)*(f_sample/N);

% Plot the frequency spectrum after filtering
figure;
plot(freqs(1:N/2), abs(frequency_content_filtered(1:N/2)));
title('FFT Analysis After Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([f_transmit-100000 f_transmit+100000]);  % Zoom in to Doppler shift region

% Generate spectrogram to visualize the cleaned Doppler signal
figure;
spectrogram(filtered_signal_final, 256, 200, 512, f_sample, 'yaxis');
title('Spectrogram of Filtered Doppler Signal');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% Mark the detected Doppler frequency
[max_amp, idx] = max(abs(frequency_content_filtered(1:N/2)));
detected_freq = freqs(idx);
hold on;
stem(detected_freq, max_amp, 'r');
%legend('FFT Spectrum', 'Doppler Shift');
fprintf('Filtered Estimated Doppler Shift Frequency: %.2f Hz\n', detected_freq);

% Compare SNR before and after denoising
snr_before = snr(received_signal_noisy);
snr_after = snr(filtered_signal_final);
fprintf('SNR Before Denoising: %.2f dB\n', snr_before);
fprintf('SNR After Denoising: %.2f dB\n', snr_after);

figure;
subplot(2,1,1);
plot(time(1:L), received_signal_noisy(1:L));
title('Received Noisy Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(time(1:L), filtered_signal_final(1:L));
title('Denoised Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Estimate velocity from detected Doppler shift
f_estimated_shift = detected_freq - f_transmit;
v_estimated = (f_estimated_shift * c) / (2 * f_transmit * cos(theta));

fprintf('Estimated Blood Flow Velocity: %.2f m/s\n', v_estimated);

angles = [30, 45, 60];  % Different angles for evaluation
v_corrected = zeros(size(angles));

for i = 1:length(angles)
    theta_i = deg2rad(angles(i));
    v_corrected(i) = (f_estimated_shift * c) / (2 * f_transmit * cos(theta_i));
    fprintf('Estimated Velocity at %.0f degrees: %.2f m/s\n', angles(i), v_corrected(i));
end

% Plot velocity estimation for different angles
figure;
plot(angles, v_corrected, '-o');
title('Velocity Estimation vs. Angle');
xlabel('Angle (degrees)');
ylabel('Estimated Velocity (m/s)');
grid on;