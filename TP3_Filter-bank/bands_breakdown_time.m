clear all; clc; close all;

%% Chargement du signal
file_path = 'pcm_48k.mat'; % Indiquez le bon chemin vers le fichier .mat
if exist(file_path, 'file') == 2
    load(file_path, '-mat'); % Chargement du fichier
else
    error('Fichier non trouvé : %s', file_path);
end

% Vérification des données
if ~exist('pcm_48', 'var')
    error('La variable pcm_48 n''existe pas dans le fichier chargé.');
end

Xin = pcm_48(5e4:1e5); % Signal d'entrée
Fs = 48000;   % Fréquence d'échantillonnage
N = length(Xin);
t = (0:N-1) / Fs;

%% Frequencies for octave bands (20 Hz to 20 kHz)
f_min = 20;                   % Minimum frequency
num_bands = 10;               % Number of bands
frequencies = f_min * 2.^(0:(num_bands-1)); % Octave band center frequencies

%% Filter design
fir_order = 2000;            % FIR filter order
filters = cell(1, num_bands); % Filter storage
delays = zeros(1, num_bands); % Delay storage
for i = 1:num_bands
    if i == 1
        % Low-pass filter for the lowest band (up to 40 Hz)
        cutoff = frequencies(i+1) / (Fs / 2); % Normalized cutoff frequency
        filters{i} = fir1(fir_order, cutoff, 'low'); % Low-pass filter
    elseif i == num_bands
        % High-pass filter for the highest band (beyond 10 kHz)
        cutoff = frequencies(i) / (Fs / 2); % Normalized cutoff frequency
        filters{i} = fir1(fir_order, cutoff, 'high'); % High-pass filter
    else
        % Band-pass filter for intermediate bands
        cutoff = [frequencies(i), frequencies(i+1)] / (Fs / 2); % Normalized cutoff
        filters{i} = fir1(fir_order, cutoff, 'bandpass'); % Band-pass filter
    end
    % Compute and store group delay (filter order / 2)
    delays(i) = fir_order / 2;
end

% Apply filters and correct for delay
subbands = cell(1, num_bands);
for i = 1:num_bands
    % Filter signal
    filtered_signal = filter(filters{i}, 1, Xin);
    % Compensate for delay by circular shift
    subbands{i} = circshift(filtered_signal, -delays(i));
end

% Reconstruct signal by summing subbands
y_reconstructed = zeros(size(Xin));
for i = 1:num_bands
    y_reconstructed = y_reconstructed + subbands{i};
end

%% Visualize signals and errors
figure;
ploted_samples = 500;
x_lim = [0 (ploted_samples)/Fs];

% (1) Original and reconstructed signal
subplot(2,1,1);
plot(t, Xin, 'b', 'LineWidth', 1.5); hold on;
plot(t, y_reconstructed, 'r', 'LineWidth', 1.5);
legend('Original Signal', 'Reconstructed Signal');
title('Original vs Reconstructed Signal');
xlabel('Time (s)');
ylabel('Amplitude');
xlim(x_lim);
grid on;

% (2) Reconstruction error
subplot(2,1,2);
error_signal = Xin - y_reconstructed;
plot(t, error_signal, 'k', 'LineWidth', 1.5);
title('Reconstruction Error');
xlabel('Time (s)');
ylabel('Amplitude');
xlim(x_lim);
ylim([-2e-4 4e-4])
grid on;

% (3) Bode diagram for all filters
figure;
hold on;
for i = 1:num_bands
    [H, f] = freqz(filters{i}, 1, 4096, Fs); % Frequency response
    semilogx(f, 20*log10(abs(H)), 'LineWidth', 1.5);
end
title('Bode Diagram of Octave Band Filters');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([10e1 1e5]);
grid on;
legend(arrayfun(@(x) sprintf('Band %d', x), 1:num_bands, 'UniformOutput', false));
hold off;

% Overall figure title
sgtitle('Analysis of Octave Bands (20 Hz to 20 kHz)');

