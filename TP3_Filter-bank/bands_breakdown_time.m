clear all; clc; close all;

%% Load Signal
Xin = load("pcm_48k.mat", '-mat').pcm_48; % Input PCM signal
Fs = 48000;                              % Sampling frequency
N = length(Xin);                         % Signal length
t = (0:N-1)/Fs;                          % Time vector

%% FIR Perfect Reconstruction Filter Bank Design
n  = 21;          % Filter order (must be odd)
fp = 0.25;         % Initial passband edge for splitting at Nyquist

% Create filters for each stage
[h0, h1, g0, g1] = firpr2chfb(n, fp); % Stage-1 filters

%% Octave Band Decomposition using 6 Bands
% Stage 1: Split signal into Lowpass and Highpass
X1_lp = filter(h0, 1, Xin); % Lowpass band
X1_hp = filter(h1, 1, Xin); % Highpass band

% Downsample both bands
X1_lp_ds = downsample(X1_lp, 2);
X1_hp_ds = downsample(X1_hp, 2);

% Stage 2: Split Lowpass (X1_lp_ds) into two subbands
X2_lp = filter(h0, 1, X1_lp_ds);
X2_hp = filter(h1, 1, X1_lp_ds);

X2_lp_ds = downsample(X2_lp, 2);
X2_hp_ds = downsample(X2_hp, 2);

% Stage 3: Split Highpass from Stage 1 into two subbands
X3_lp = filter(h0, 1, X1_hp_ds);
X3_hp = filter(h1, 1, X1_hp_ds);

X3_lp_ds = downsample(X3_lp, 2);
X3_hp_ds = downsample(X3_hp, 2);

%% Store all 6 Bands (already downsampled)
Band1 = X2_lp_ds; % Lowest band
Band2 = X2_hp_ds; % Low-mid band
Band3 = X3_lp_ds; % Mid band
Band4 = X3_hp_ds; % High-mid band
Band5 = downsample(filter(h0, 1, X3_hp_ds), 2); % Subband splitting continues
Band6 = downsample(filter(h1, 1, X3_hp_ds), 2); % Highest band

%% Reconstruct Signal
% Upsample and reconstruct at each stage
Y3_lp_us = upsample(Band3, 2);
Y3_hp_us = upsample(Band4, 2);
Y3_lp = filter(g0, 1, Y3_lp_us);
Y3_hp = filter(g1, 1, Y3_hp_us);
Y3 = Y3_lp + Y3_hp;

Y2_lp_us = upsample(Band1, 2);
Y2_hp_us = upsample(Band2, 2);
Y2_lp = filter(g0, 1, Y2_lp_us);
Y2_hp = filter(g1, 1, Y2_hp_us);
Y2 = Y2_lp + Y2_hp;

Y_final_lp_us = upsample(Y2, 2);
Y_final_hp_us = upsample(Y3, 2);
Y_final_lp = filter(g0, 1, Y_final_lp_us);
Y_final_hp = filter(g1, 1, Y_final_hp_us);

X_reconstructed = Y_final_lp + Y_final_hp;

%% Frequency Response of Octave Bands
figure;
N_fft = 2048;
f = linspace(0, Fs/2, N_fft); % Frequency axis

% Plot FFT of each band
hold on;
plot(f, abs(fft(Band1, N_fft)), 'LineWidth', 1.5, 'DisplayName', 'Band 1');
plot(f, abs(fft(Band2, N_fft)), 'LineWidth', 1.5, 'DisplayName', 'Band 2');
plot(f, abs(fft(Band3, N_fft)), 'LineWidth', 1.5, 'DisplayName', 'Band 3');
plot(f, abs(fft(Band4, N_fft)), 'LineWidth', 1.5, 'DisplayName', 'Band 4');
plot(f, abs(fft(Band5, N_fft)), 'LineWidth', 1.5, 'DisplayName', 'Band 5');
plot(f, abs(fft(Band6, N_fft)), 'LineWidth', 1.5, 'DisplayName', 'Band 6');

title('Octave Bands in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('show');
grid on;

%% Compare Original and Reconstructed Signal
figure;

subplot(3,1,1);
plot(t, Xin, 'b');
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(t, X_reconstructed, 'r');
title('Reconstructed Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,3);
plot(t, Xin - X_reconstructed, 'k');
title('Reconstruction Error');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% Display Reconstruction Error
error = max(abs(Xin - X_reconstructed));
disp(['Maximum Reconstruction Error: ', num2str(error)]);
