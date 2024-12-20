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

%% First Cascade Decimation
% Decimation Parameters
Fs_pcm = 48e3;  % Target sampling frequency
M1 = 2;         % First decimation factor
M2 = 8;         % Second decimation factor
M3 = 8;         % Third decimation factor

%% Cascade Decimation
pdm_decimated = decimate_cascade(pdm_in, Fs_pdm, M1, M2, M3);

%% FFT of Decimated Signal
PDM_dec_fft = abs(fft(pdm_decimated)); % FFT Magnitude
PDM_dec_fft = PDM_dec_fft / max(PDM_dec_fft); % Normalization
PDM_dec_fft_db = 20 * log10(PDM_dec_fft); % Convert to dB

N_dec = length(PDM_dec_fft); % Number of samples
freq_dec = (0:N_dec-1)*(Fs_pcm/N_dec); % Frequency axis for the decimated signal

figure(2)
subplot(2,1,1)
plot(freq_dec/1e6, PDM_dec_fft_db)
xlim([0 Fstop/1e6])
grid on
hold on

subplot(2,1,2)
plot(pdm_decimated(1:300))
grid on
hold on

%% Second cascade Decimation
% Decimation Parameters
Fs_pcm = 48e3;  % Target sampling frequency
M1 = 2;         % First decimation factor
M2 = 4;         % Second decimation factor
M3 = 4;         % Third decimation factor
M4 = 4;         % Fourth decimation factor

%% Cascade Decimation
pdm_decimated = decimate_cascade(pdm_in, Fs_pdm, M1, M2, M3, M4);

%% FFT of Decimated Signal
PDM_dec_fft = abs(fft(pdm_decimated)); % FFT Magnitude
PDM_dec_fft = PDM_dec_fft / max(PDM_dec_fft); % Normalization
PDM_dec_fft_db = 20 * log10(PDM_dec_fft); % Convert to dB

N_dec = length(PDM_dec_fft); % Number of samples
freq_dec = (0:N_dec-1)*(Fs_pcm/N_dec); % Frequency axis for the decimated signal

subplot(2,1,1)
plot(freq_dec/1e6, PDM_dec_fft_db)
xlim([0 Fstop/1e6])
grid on

subplot(2,1,2)
plot(pdm_decimated(1:300))
grid on

%% Third cascade Decimation
% Decimation Parameters
Fs_pcm = 48e3;  % Target sampling frequency
M1 = 8;         % First decimation factor
M2 = 4;         % Second decimation factor
M3 = 4;         % Third decimation factor

%% Cascade Decimation
pdm_decimated = decimate_cascade(pdm_in, Fs_pdm, M1, M2, M3);

%% FFT of Decimated Signal
PDM_dec_fft = abs(fft(pdm_decimated)); % FFT Magnitude
PDM_dec_fft = PDM_dec_fft / max(PDM_dec_fft); % Normalization
PDM_dec_fft_db = 20 * log10(PDM_dec_fft); % Convert to dB

N_dec = length(PDM_dec_fft); % Number of samples
freq_dec = (0:N_dec-1)*(Fs_pcm/N_dec); % Frequency axis for the decimated signal

subplot(2,1,1)
plot(freq_dec/1e6, PDM_dec_fft_db)
xlim([0 Fstop/1e6])
grid on

subplot(2,1,2)
plot(pdm_decimated(1:300))
grid on

%% Fourth cascade Decimation
% Decimation Parameters
Fs_pcm = 48e3;  % Target sampling frequency
M1 = 4;         % First decimation factor
M2 = 4;         % Second decimation factor
M3 = 8;         % Third decimation factor
M4 = 4;         % Fourth decimation factor
M5 = 2;         % Fourth decimation factor

%% Cascade Decimation
pdm_decimated = decimate_cascade(pdm_in, Fs_pdm, M1, M2, M3);

%% FFT of Decimated Signal
PDM_dec_fft = abs(fft(pdm_decimated)); % FFT Magnitude
PDM_dec_fft = PDM_dec_fft / max(PDM_dec_fft); % Normalization
PDM_dec_fft_db = 20 * log10(PDM_dec_fft); % Convert to dB

N_dec = length(PDM_dec_fft); % Number of samples
freq_dec = (0:N_dec-1)*(Fs_pcm/N_dec); % Frequency axis for the decimated signal

subplot(2,1,1)
plot(freq_dec/1e6, PDM_dec_fft_db)
xlabel('Frequency (MHz)')
ylabel('Normalized Amplitude (dB)')
title('FFT of Cascade Decimations')
legend('2x8x8', '2x4x4x4', '8x4x4', '4x4x8')
xlim([0 Fstop/1e6])
grid on

subplot(2,1,2)
plot(pdm_decimated(1:300))
xlabel('Sample Index')
ylabel('Amplitude')
title('Cascade Decimated PCM')
legend('2x8x8', '2x4x4x4', '8x4x4', '4x4x8')
axis([0 100 -1.1 1.1])
grid on

%% Save the Decimated Output
% audiowrite("pcm_from_pdm.wav", pdm_decimated, Fs_pcm);