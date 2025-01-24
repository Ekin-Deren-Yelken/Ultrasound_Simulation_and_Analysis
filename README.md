# Ultrasound_Simulation_and_Analysis
Three projects in teh simulation and analysis of python

Ultrasound Simulation and Analysis
Projects:

1. Write a simulation of ultrasound wave propagation using finite element analysis (FEA).
2. Implement an algorithm to calculate Doppler shift for blood flow measurements.
3. Develop a software tool to visualize ultrasound transducer array patterns.
Ideal Language: Python (for simulation and visualization), MATLAB (for signal modeling), C++ (for performance optimization)


# Doppler Shift Calculation for Blood Flow Measurements

## Overview
This project focuses on developing an algorithm to calculate the Doppler shift in ultrasound signals based on interactions with moving blood particles, enabling blood flow velocity estimation.

### Objectives
- Simulate ultrasound wave interactions with moving particles.
- Implement signal processing techniques to extract Doppler frequency shifts.
- Visualize and analyze the results using  MATLAB.

---

## Background

### Doppler Effect in Ultrasound Imaging
The Doppler effect describes the frequency shift that occurs when ultrasound waves reflect off moving blood cells. The relationship between the transmitted frequency, received frequency, and velocity of moving particles is given by the following equation:

$f_d = \frac{2 f_t v \cos(\theta)}{c}$

Where:  
- $f_d$ = Doppler shift (Hz)  
- $f_t$ = Transmitted frequency (Hz)  
- $v$ = Blood flow velocity (m/s)  
- $&theta;$ = Angle between ultrasound beam and blood flow direction (degrees)  
- $c$ = Speed of sound in biological tissue (~1540 m/s)  

### Types of Doppler Ultrasound
- **Continuous Wave (CW) Doppler:** Constant wave transmission and reception for high-velocity detection.
- **Pulsed Wave (PW) Doppler:** Short bursts of ultrasound to measure velocity at specific depths.

---

### Transmitted Signal
This figure shows the clean transmitted ultrasound signal at $f = 2 MHz$. The signal represents the ideal waveform generated by the ultrasound transducer.

#### Key Observations:
- Pure sinusoidal waveform at the transmitted frequency.
- No external noise or Doppler shift applied.

$1e-03$ second signal

![transmitted-signal](https://github.com/user-attachments/assets/85c70544-94c3-44dd-9a03-26b5eea79410)


### Noise
The received signal includes the effects of Doppler shift, Gaussian noise, powerline interference (50 Hz), and motion artifacts (0.5 Hz).

#### Key Observations:
- Amplitude modulations due to pulsatile blood flow.
- High-frequency noise and low-frequency motion artifacts present.

![recieved-signal](https://github.com/user-attachments/assets/d3503c4c-91c4-48c9-9384-d013ee14f7a4)

### Frequecy Spectrum of Noisy Signal
The FFT analysis of the noisy received signal showing prominent frequency components, including the Doppler shift and unwanted noise at 50 Hz.

#### Key Observations:
- A clear peak near 2 MHz corresponding to the Doppler frequency.
- Additional peaks representing noise components such as powerline interference, namely close to 0.1 MHz and 4 MHz

![FFT-recieved-w-doppler](https://github.com/user-attachments/assets/13b64573-44b9-4b2c-949a-51f66da4d795)

Calcualted Doppler Shift Frequency (Mean): 919.76 Hz

### Spectrogram of Noisy Signal
A time-frequency representation of the received signal, illustrating how different frequency components evolve over time.

#### Key Observations:
- Doppler frequency shift is visible over time.
- Noise components persist throughout the signal duration.

![spectrogram-noisy](https://github.com/user-attachments/assets/c8129ba2-d032-4ba1-8efa-a5101cdc1e0f)

### Filtering
The received signal after applying notch filtering (to remove 50 Hz powerline noise), high-pass filtering (to remove motion artifacts), and low-pass filtering (to suppress high-frequency noise).

![recieved-filter](https://github.com/user-attachments/assets/e24b0cf1-84e3-439b-990e-de04f0130dad)

### Frequency Spectrum and Spectrogram of Filtered Signal
FFT analysis of the filtered signal, showing suppression of unwanted noise while retaining the Doppler frequency component.

####Key Observations:
- The 50 Hz powerline noise and high-frequency components are effectively removed.
- The Doppler shift peak is the only remaining distinguishable peak.
- Spectrogram is significantly clearer.

![fft-recieved-after-filter](https://github.com/user-attachments/assets/bcee7f88-a887-4f99-b5c8-ec240c4788f3)

![spectrogram-after-filter](https://github.com/user-attachments/assets/0120abf7-3860-4090-8c10-c8c223868507)

### Denoising Efficacy
After denoising, the doppler shift frequency is estimated at 2001000 Hz or 2.001 MHz. This error can be attributatble to memory errors in operating with large numbers on computers. Comparing the two signals, we can calculate the signal to noise ration before and after filtering to get a beter idea of how well teh filters worked.

![filtered-compare](https://github.com/user-attachments/assets/3fd5e53d-73b4-4048-98dc-2669c774adc2)

SNR before filtering = 0.28 dB
SNR after filtering = 2.09 dB

The blood flow is estimated after denoising and tested at different angles:
![image](https://github.com/user-attachments/assets/0363a12a-900b-4a42-a612-b1dc1b2ea3c8)
