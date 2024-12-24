clear all; clc; close all;

% Chargement du signal
load('pcm_48k.mat', '-mat'); % Chargement du signal d'entrée
Xin = pcm_48;  % Signal d'entrée
Fs = 48000;    % Fréquence d'échantillonnage
N = length(Xin);
t = (0:N-1)/Fs;

% Définir les bandes d'octaves
f_min = 20; % Fréquence minimale (en Hz)
f_max = 20000; % Fréquence maximale (Nyquist)
num_bands = 10;
center_freqs = logspace(log10(f_min), log10(f_max), num_bands);

% Design des filtres FIR
filters = cell(1, num_bands);
for i = 1:num_bands
    if i == 1
        % Filtre passe-bas pour la première bande
        f_high = center_freqs(1);
        f_high_norm = f_high / (Fs/2);

        % FIR passe-bas avec fenêtre de Hamming
        n = 200; % Ordre du filtre (fixe pour simplicité)
        filters{i} = fir1(n, f_high_norm, 'low', hamming(n+1));
    elseif i == num_bands
        % Filtre passe-haut pour la dernière bande
        f_low = center_freqs(num_bands-1);
        f_low_norm = f_low / (Fs/2);

        % FIR passe-haut avec fenêtre de Hamming
        n = 200; % Ordre du filtre
        filters{i} = fir1(n, f_low_norm, 'high', hamming(n+1));
    else
        % Filtre passe-bande pour les bandes intermédiaires
        f_low = center_freqs(i-1);
        f_high = center_freqs(i);
        f_low_norm = f_low / (Fs/2);
        f_high_norm = f_high / (Fs/2);

        % FIR passe-bande avec fenêtre de Hamming
        n = 200; % Ordre du filtre
        filters{i} = fir1(n, [f_low_norm f_high_norm], 'bandpass', hamming(n+1));
    end
end

% Décomposer le signal en sous-bandes
subbands = zeros(num_bands, N);
for i = 1:num_bands
    b = filters{i};
    subbands(i, :) = filter(b, 1, Xin); % Application du filtre FIR
end

% Reconstruction du signal
reconstructed_signal = sum(subbands, 1);

% Gestion des retards des filtres FIR (délais linéaires)
group_delays = zeros(1, num_bands);
for i = 1:num_bands
    group_delays(i) = length(filters{i}) / 2; % Délais linéaires pour filtres FIR
end
max_delay = max(group_delays);
reconstructed_signal = circshift(reconstructed_signal, -round(max_delay)); % Compense le retard

% Calcul de l'erreur
error_signal = Xin - reconstructed_signal;

% Visualisation
figure;
subplot(3, 1, 1);
plot(t, Xin);
title("Signal d'origine");
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(t, reconstructed_signal);
title("Signal reconstruit");
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(3, 1, 3);
plot(t, error_signal);
title("Erreur entre le signal original et reconstruit");
xlabel("Temps (s)");
ylabel("Amplitude");

% Affichage des réponses fréquentielles des filtres
figure;
hold on;
for i = 1:num_bands
    b = filters{i};
    [H, f] = freqz(b, 1, 1024, Fs);
    plot(f, 20*log10(abs(H)));
end
title("Réponses fréquentielles des filtres FIR");
xlabel("Fréquence (Hz)");
ylabel("Amplitude (dB)");
grid on;
legend(arrayfun(@(f) sprintf("Bande %d", f), 1:num_bands, "UniformOutput", false));

