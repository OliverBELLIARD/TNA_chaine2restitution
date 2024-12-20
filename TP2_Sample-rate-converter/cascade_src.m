%% Sample-Rate Conversion: 44.1kHz to 48kHz (Cascade Interpolation)
% Conversion utilisant une cascade d'interpolation/decimation
clear all; clc; close all;

%% Parameters for Sample Rate Converter
Fs_in  = 44100;   % Fréquence d'entrée (Hz)
Fs_out = 48000;   % Fréquence de sortie (Hz)

% Chargement du signal d'entrée
input_signal = load('playback_44100.mat', "-mat").w441(8e4:1e5);
N = length(input_signal);
t_in = (0:N-1) / Fs_in; % Vecteur temps du signal d'entrée

% Interpolation/Decimation facteurs en cascade
cascade_factors = [
    8, 7;  % Étape 1: 8/7
    5, 3;  % Étape 2: 5/3
    4, 7   % Étape 3: 4/7
];
num_stages = size(cascade_factors, 1); % Nombre d'étapes

%% Étape par Étape: Interpolation et Décimation
intermediate_signals = cell(num_stages + 1, 1);
intermediate_signals{1} = input_signal; % Stockage du signal d'entrée initial

current_signal = input_signal; % Signal courant
current_fs = Fs_in; % Fréquence actuelle

for stage = 1:num_stages
    L = cascade_factors(stage, 1); % Facteur d'interpolation
    M = cascade_factors(stage, 2); % Facteur de décimation
    fprintf('Stage %d: Interpolation %d / Decimation %d\n', stage, L, M);

    % Interpolation
    interpolated_signal = interpolate_audio(current_signal, current_fs, L, 16);

    % Décimation
    current_fs = (current_fs * L); % Nouvelle fréquence d'échantillonnage
    current_signal = decimate_audio(interpolated_signal, current_fs, M, 16);

    current_fs = (current_fs) / M; % Nouvelle fréquence d'échantillonnage
    % Stockage du signal intermédiaire
    intermediate_signals{stage + 1} = current_signal;
end

%% Signal final
out_signal = current_signal;
t_out = 0:1/Fs_out:(length(out_signal)-1)/Fs_out; % Temps final

%% Comparaison Étape par Étape
figure;
for stage = 1:num_stages
    subplot(2, 1, 1);
    title('Signaux intermédiaires');
    hold on;
    t = (0:length(intermediate_signals{stage})-1) / Fs_in;
    plot(t, intermediate_signals{stage});
    xlabel('Temps (s)');
    ylabel('Amplitude');
    legend('Etage 8/7', 'Etage 5/3', 'Etage 4/7')
    grid on;

    subplot(2, 1, 2);
    hold on;
    N_fft = length(intermediate_signals{stage});
    frequencies = linspace(0, current_fs/2, N_fft/2);
    fft_signal = abs(fft(intermediate_signals{stage}, N_fft));
    semilogx(frequencies, 20 * log10(fft_signal(1:N_fft/2)));
    title('Réponses Fréquentielles intermédiaires');
    xlabel('Fréquence (Hz)');
    ylabel('Amplitude (dB)');
    legend('Etage 8/7', 'Etage 5/3', 'Etage 4/7')
    grid on;
end
hold off;

%% Comparaison Signal d'Entrée vs Signal Final
figure;
subplot(2, 1, 1);
plot(t_in, input_signal, 'b', 'DisplayName', 'Signal Entrée (44.1 kHz)');
hold on;
plot(t_out, out_signal, 'r', 'DisplayName', 'Signal Final (48 kHz)');
hold off;
title("Comparaison: Signal d\'Entrée et Signal Final");
xlabel('Temps [s]');
ylabel('Amplitude');
legend;
grid on;

subplot(2, 1, 2);
N_fft = max(length(input_signal), length(out_signal));
frequencies_in = linspace(0, Fs_in/2, N_fft/2);
frequencies_out = linspace(0, Fs_out/2, N_fft/2);

input_fft = abs(fft(input_signal, N_fft));
output_fft = abs(fft(out_signal, N_fft));

plot(frequencies_in, 20 * log10(input_fft(1:N_fft/2)), 'b', 'DisplayName', 'Signal Entrée FFT');
hold on;
plot(frequencies_out, 20 * log10(output_fft(1:N_fft/2)), 'r', 'DisplayName', 'Signal Final FFT');
hold off;
title('Comparaison: Réponse Fréquentielle');
xlabel('Fréquence (Hz)');
ylabel('Amplitude [dB]');
legend;
grid on;

%% Sauvegarde des Fichiers Audio
audiowrite('input_signal_44k.wav', input_signal, Fs_in);
audiowrite('output_signal_48k.wav', out_signal, Fs_out);

% Optionnel: Jouer le signal
% sound(input_signal, Fs_in);
% sound(out_signal, Fs_out);