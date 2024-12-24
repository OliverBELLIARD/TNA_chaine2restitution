clear all; clc; close all;

% Chargement du signal
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

% Définir les bandes d'octaves
f_min = 20;    % Fréquence minimale (Hz)
f_max = 20000; % Fréquence maximale (Hz)
num_bands = 10;
center_freqs = logspace(log10(f_min), log10(f_max), num_bands);

% Design des filtres passe-bande
filters = cell(1, num_bands);
for i = 1:num_bands
    if i == 1
        % Passe-bas pour la première bande
        f_high = center_freqs(1);
        f_high_norm = f_high / (Fs / 2);
        n = 200; % Ordre du filtre
        filters{i} = fir1(n, f_high_norm, 'low');
    elseif i == num_bands
        % Passe-haut pour la dernière bande
        f_low = center_freqs(num_bands - 1);
        f_low_norm = f_low / (Fs / 2);
        n = 200; % Ordre du filtre
        filters{i} = fir1(n, f_low_norm, 'high');
    else
        % Passe-bande pour les bandes intermédiaires
        f_low = center_freqs(i - 1);
        f_high = center_freqs(i);
        f_low_norm = f_low / (Fs / 2);
        f_high_norm = f_high / (Fs / 2);
        n = 200; % Ordre du filtre
        filters{i} = fir1(n, [f_low_norm f_high_norm], 'bandpass');
    end
end

% Application des filtres
subbands = zeros(num_bands, N);
for i = 1:num_bands
    subbands(i, :) = filter(filters{i}, 1, Xin);
end

% Reconstruction du signal
reconstructed_signal = sum(subbands, 1);

% Calcul de l'erreur
error_signal = Xin - reconstructed_signal;

% Visualisation
figure;
subplot(3, 1, 1);
plot(t, Xin);
title('Signal d''origine');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(t, reconstructed_signal);
title('Signal reconstruit');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(3, 1, 3);
plot(t, error_signal);
title('Erreur');
xlabel('Temps (s)');
ylabel('Amplitude');

