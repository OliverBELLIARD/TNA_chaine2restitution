%% Sample-Rate Conversion: 44.1kHz to 48kHz
% This script designs and implements a Sample-Rate Converter (SRC)
% to convert a signal sampled at 44.1 kHz to 48 kHz.

clear all; clc; close all;

%% Parameters for the Sample Rate Converter
Fs_in  = 44100;  % Input sampling rate (Hz)
Fs_out = 48000;  % Target sampling rate (Hz)

% L and M are the interpolation and decimation factors
[L, M] = rat(Fs_out / Fs_in); % Rational approximation for resampling ratio

%% Load the input signal
input_signal = load('playback_44100.mat', "-mat").w441(50e3:1e5);
N = length(input_signal);
t_in = (0:N-1)/Fs_in; % Time vector for input signal

%% Interpolation (Upsampling)
% Anti-Aliasing Filter Design
Nbits = 16;
Fp = 20e3;           % Passband edge frequency
Fst = Fs_in / 2;     % Stopband edge frequency
Fs_up = Fs_in * L;
Apass = 0.01;        % Passband ripple in dB
Astop = 1.76 + 6.02 * 16; % Stopband attenuation in dB

interpolated_signal = interpolate_audio(input_signal, Fs_in, L, Nbits);

%% Decimation (Downsampling)
out_signal = decimate_audio(interpolated_signal, Fs_in * L, M, Nbits);

%% Time Vectors for Outputs
t_out = 0:1/Fs_out:(length(out_signal)-1)/Fs_out;

%% Plot the Input Signal
figure;
subplot(2,1,1);
hold on
plot(t_in, input_signal);

%% Plot the Converted Signal
plot(t_out, out_signal);
title('Input (44.1 kHz) and Output (48 kHz) Signals');
xlabel('Time [s]'); ylabel('Amplitude'); grid on;
hold off

%% Verify Frequency Content via FFT
N_fft = N;
frequencies_in  = linspace(0, Fs_in/2, N_fft/2);
frequencies_out = linspace(0, Fs_out/2, N_fft/2);

input_fft = abs(fft(input_signal, N_fft));
output_fft = abs(fft(out_signal, N_fft));

input_fft_db = 20 * log10(input_fft);
output_fft_db = 20 * log10(output_fft);

figure;
subplot(2,1,2);
hold on
loglog(frequencies_in, 20*log10(input_fft(1:N_fft/2)));
loglog(frequencies_out, 20*log10(output_fft(1:N_fft/2)));
hold off
title('Frequency Response of Input (44.1 kHz) and output Signal (48 kHz)');
xlabel('Frequency [Hz]'); ylabel('Amplitude [dB]'); grid on;
legend('Input (44.1 kHz)', 'output (48 kHz)');

%% Audio Output
% Optionally save or play the audio files
audiowrite('input_signal_44k.wav', input_signal, Fs_in);
audiowrite('output_signal_48k.wav', out_signal, Fs_out);
sound(input_signal, Fs_in);
sound(out_signal, Fs_out);