clear all; clc; close all;

%% Parameters for the Sample Rate Converter
Fs_in  = 44100;  % Input sampling rate (Hz)
Fs_out = 48000;  % Target sampling rate (Hz)

% L and M are the interpolation and decimation factors
[L, M] = rat(Fs_out / Fs_in); % Rational approximation for resampling ratio

%% Input signal
input_signal = load('playback_44100.mat', "-mat").w441(50e3:1e5);
% Audio file creation
% audiowrite("playback_44100.wav", input_signal, Fs_in);
% sound(input_signal, Fs_in);
N = length(input_signal);

%% Plot the Input Signal
figure;
subplot(2,1,1);
plot(input_signal);
title('Input Signal (44.1 kHz)');
xlabel('Time [s]'); ylabel('Amplitude'); grid on;

%% Design the Anti-Aliasing Filter
Fp = min(Fs_in, Fs_out) / 2; % Passband edge frequency
Fst = Fp + 1000; % Stopband edge frequency
Apass = 0.01; % Passband ripple in dB
Astop = 80;   % Stopband attenuation in dB

h = fdesign.lowpass('Fp,Fst,Ap,Ast', Fp, Fst, Apass, Astop, Fs_in * L);
anti_alias_filter = design(h, 'ellip', 'SystemObject', true);

%% Perform the Interpolation (Upsampling + Filtering)
upsampled_signal = zeros(1, L * length(input_signal));
upsampled_signal(1:L:end) = input_signal; % Upsampling
filtered_signal = step(anti_alias_filter, upsampled_signal);

%% Decimation (Downsampling)
out_signal = filtered_signal(1:M:end);

%% Time Vectors for Outputs
t_out = 0:1/Fs_out:(length(out_signal)-1)/Fs_out;

%% Audio Output
% Optionally save or play the audio files
audiowrite('input_signal_44k.wav', input_signal, Fs_in);
audiowrite('output_signal_48k.wav', out_signal, Fs_out);
% sound(out_signal, Fs_out);

%% Plot the Converted Signal
subplot(2,1,2);
plot(t_out, out_signal);
title('Output Signal (48 kHz)');
xlabel('Time [s]'); ylabel('Amplitude'); grid on;

%% Verify Frequency Content via FFT
N_fft = 4096;
frequencies_in  = linspace(0, Fs_in/2, N_fft/2);
frequencies_out = linspace(0, Fs_out/2, N_fft/2);

input_fft = abs(fft(input_signal, N_fft));
output_fft = abs(fft(out_signal, N_fft));

figure;
subplot(2,1,1);
plot(frequencies_in, 20*log10(input_fft(1:N_fft/2)));
title('Frequency Response of Input Signal (44.1 kHz)');
xlabel('Frequency [Hz]'); ylabel('Amplitude [dB]'); grid on;

subplot(2,1,2);
plot(frequencies_out, 20*log10(output_fft(1:N_fft/2)));
title('Frequency Response of Output Signal (48 kHz)');
xlabel('Frequency [Hz]'); ylabel('Amplitude [dB]'); grid on;
