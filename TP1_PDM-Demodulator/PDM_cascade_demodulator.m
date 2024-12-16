clear all; clc; close all;

pdm_in = load('pdm_in.mat', "-mat").in;
Fs_pdm = 6.144e6;
N = length(pdm_in); % Number of samples

% Frequency axis
freq = (0:N-1)*(Fs_pdm/N); % Frequency axis for single-sided FFT

%% Son du signal PDM d'entrée
% Création du fichier audio
audiowrite("pdm_in.wav", pdm_in, Fs_pdm);

%% FFT of the PDM input signal
PDM_fft = abs(fft(pdm_in)); % Magnitude of the FFT
PDM_fft = PDM_fft / max(PDM_fft); % Normalization

% Convert to dB
PDM_fft_db = 20 * log10(PDM_fft); % Convert amplitude to decibels (log scale)

%% Visualisation du signal PDM d'entrée et FFT
figure(1)

% FFT plot in dB (frequency domain)
subplot(2,1,1)
plot(freq/1e6, PDM_fft_db) % Frequency in MHz
xlabel('Frequency (MHz)')
ylabel('Normalized Amplitude')
title('FFT of PDM Input Signal')
%xlim([0 Fs_pdm/2e6]) % Plot up to Nyquist frequency
xlim([0 0.03]) % Plot up to noise frequencies
grid on

% PDM signal (time domain)
subplot(2,1,2)
plot(pdm_in(1:300))
xlabel('Sample Index')
ylabel('Amplitude')
title('PDM Input Signal (Time Domain)')
axis([0 300 -1.1 1.1])
grid on

%% Filter desing
Fs_pcm = 48e3;
M = floor(Fs_pdm/Fs_pcm); % Decimation factor
 
Apass = power(10, 0.1/20); % Filter gain
Fpass = 20e3; % Cut off frequency
Fstop = Fs_pdm/(2*M); % Frequency of the attenuated band
Astop = 1.76 + 6.02*20; % Gain in the attenuated band

%% Filtering and Decimation
% Apply cascade decimation
M1 = 16;
M2 = 4;
M3 = 2;

Fs1 = floor(Fs_pdm/M1); % New decimation factor
pdm_filtered1 = decimate_audio(pdm_in, Fs_pdm, M1);

Fs2 = floor(Fs_pdm/M2); % New decimation factor
pdm_filtered2 = decimate_audio(pdm_filtered1, Fs1, M2);

Fs3 = floor(Fs_pdm/M3); % New decimation factor
pdm_decimated = decimate_audio(pdm_filtered2, Fs2, M3);

%% FFT of the decimated input signal
PDM_dec_fft = abs(fft(pdm_decimated)); % Magnitude of the FFT
PDM_dec_fft = PDM_dec_fft / max(PDM_dec_fft); % Normalization

% Convert to dB
PDM_dec_fft_db = 20 * log10(PDM_dec_fft); % Convert amplitude to decibels (log scale)

N = length(PDM_dec_fft); % Number of samples

% Frequency axis
freq = (0:N-1)*(Fs_pcm/N); % Frequency axis for single-sided FFT

%% Visualisation du signal PDM décimé et sa FFT
figure(2)

% FFT plot in dB (frequency domain)
subplot(2,1,1)
plot(freq/1e6, PDM_dec_fft_db) % Frequency in MHz
xlabel('Frequency (MHz)')
ylabel('Normalized Amplitude')
title('FFT of decimated PCM Input Signal')
xlim([0 .03]) % Plot up to Nyquist frequency
grid on

% PDM signal (time domain)
subplot(2,1,2)
plot(pdm_decimated(1:300))
xlabel('Sample Index')
ylabel('Amplitude')
title('Decimated PCM Input Signal (Time Domain)')
axis([0 300 -1.1 1.1])
grid on

audiowrite("pcm_from_pdm.wav", pdm_decimated, Fs_pcm)