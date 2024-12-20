%% Sample-Rate Conversion: 44.1kHz to 48kHz (Reordered Cascade Interpolation)
% Conversion utilisant une cascade d'interpolation/decimation avec un nouvel ordre
clear all; clc; close all;

%% Parameters for Sample Rate Converter
Fs_in  = 44100;   % Fréquence d'entrée (Hz)
Fs_out = 48000;   % Fréquence de sortie (Hz)

% Chargement du signal d'entrée
input_signal = load('playback_44100.mat', "-mat").w441(8e4:1e5);
N = length(input_signal);
t_in = (0:N-1) / Fs_in; % Vecteur temps du signal d'entrée

% Interpolation et Décimation en cascade : nouvel ordre
cascade_factors_interp = [8, 5, 4]; % Facteurs d'interpolation
cascade_factors_deci = [3, 7, 7];   % Facteurs de décimation

num_stages = length(cascade_factors_interp); % Nombre d'étapes

%% Étape par Étape: Interpolation et Décimation
intermediate_signals_interp = cell(num_stages + 1, 1);
intermediate_signals_interp{1} = input_signal; % Signal initial

current_signal = input_signal; % Signal courant
current_fs = Fs_in; % Fréquence d'échantillonnage actuelle

% Étapes d'interpolation
for stage = 1:num_stages
    L = cascade_factors_interp(stage); % Facteur d'interpolation
    fprintf('Étape %d: Interpolation %d\n', stage, L);

    % Interpolation
    interpolated_signal = interpolate_audio(current_signal, current_fs, L, 16);
    current_fs = current_fs * L; % Mise à jour de la fréquence d'échantillonnage
    current_signal = interpolated_signal;

    % Stocker le signal interpolé
    intermediate_signals_interp{stage + 1} = current_signal;
end

% Signal après interpolation complète
signal_after_interp = current_signal;
t_after_interp = (0:length(signal_after_interp)-1) / current_fs; % Vecteur temps interpolé

% Étapes de décimation
intermediate_signals_deci = cell(num_stages + 1, 1);
intermediate_signals_deci{1} = signal_after_interp;

for stage = 1:num_stages
    M = cascade_factors_deci(stage); % Facteur de décimation
    fprintf('Étape %d: Décimation %d\n', stage, M);

    % Décimation
    decimated_signal = decimate_audio(current_signal, current_fs, M, 16);
    current_fs = current_fs / M; % Mise à jour de la fréquence d'échantillonnage
    current_signal = decimated_signal;

    % Stocker le signal décimé
    intermediate_signals_deci{stage + 1} = current_signal;
end

% Signal final après décimation
out_signal = current_signal;
t_out = 0:1/Fs_out:(length(out_signal)-1)/Fs_out; % Vecteur temps final

%% Comparaison Signal d'Entrée vs Interpolé
figure;
subplot(2, 1, 1);
plot(t_in, input_signal, 'b', 'DisplayName', "Signal d'Entrée (44.1 kHz)");
hold on;
plot(t_after_interp, signal_after_interp, 'r', 'DisplayName', 'Signal Interpolé');
hold off;
title("Comparaison: Signal d'Entrée vs Signal Interpolé");
xlabel('Temps (s)');
ylabel('Amplitude');
legend;
grid on;

subplot(2, 1, 2);
N_fft = max(length(input_signal), length(signal_after_interp));
frequencies_in = linspace(0, Fs_in/2, N_fft/2);
frequencies_interp = linspace(0, (Fs_in * prod(cascade_factors_interp))/2, N_fft/2);

input_fft = abs(fft(input_signal, N_fft));
interp_fft = abs(fft(signal_after_interp, N_fft));

plot(frequencies_in, 20 * log10(input_fft(1:N_fft/2)), 'b', 'DisplayName', "Signal d'Entrée FFT");
hold on;
plot(frequencies_interp, 20 * log10(interp_fft(1:N_fft/2)), 'r', 'DisplayName', "Signal Interpolé FFT");
hold off;
title("Réponse Fréquentielle: Signal d'Entrée vs Signal Interpolé");
xlabel('Fréquence (Hz)');
ylabel('Amplitude (dB)');
legend;
grid on;

%% Comparaison Signal d'Entrée vs Signal Final
figure;
subplot(2, 1, 1);
plot(t_in, input_signal, 'b', 'DisplayName', "Signal d'Entrée (44.1 kHz)");
hold on;
plot(t_out, out_signal, 'r', 'DisplayName', 'Signal Final');
hold off;
title("Comparaison: Signal d'Entrée vs Signal Final");
xlabel('Temps (s)');
ylabel('Amplitude');
legend;
grid on;

subplot(2, 1, 2);
frequencies_out = linspace(0, Fs_out/2, N_fft/2);

output_fft = abs(fft(out_signal, N_fft));

plot(frequencies_in, 20 * log10(input_fft(1:N_fft/2)), 'b', 'DisplayName', "Signal d'Entrée FFT");
hold on;
plot(frequencies_out, 20 * log10(output_fft(1:N_fft/2)), 'r', 'DisplayName', "Signal Final FFT");
hold off;
title("Réponse Fréquentielle: Signal d'Entrée vs Signal Final");
xlabel('Fréquence (Hz)');
ylabel('Amplitude (dB)');
legend;
grid on;

%% Sauvegarde des Fichiers Audio
audiowrite('input_signal_44k.wav', input_signal, Fs_in);
audiowrite('output_signal_48k.wav', out_signal, Fs_out);

% Optionnel : Jouer les signaux
% sound(input_signal, Fs_in);
% sound(out_signal, Fs_out);
