clear; clc; close all;

% Charger le signal d'entrée
data = load('playback_44100.mat', '-mat');
input_signal = data.w441(50e3:1e5); % Utiliser une fraction du signal pour accélérer les tests
Fs_in = 44100; % Fréquence d'échantillonnage initiale
Fs_out = 48000; % Fréquence d'échantillonnage cible

%% 1. Version naïve de validation
% Calculer le PPCM des deux fréquences pour trouver le facteur de conversion exact
L = 480; % Sur-échantillonnage (multiplicateur)
M = 441; % Sous-échantillonnage (diviseur)

% Sur-échantillonnage (insertion de zéros)
signal_upsampled = upsample(input_signal, L);

% Conception d'un filtre anti-repliement simple (Filtre passe-bas à 20 kHz)
Fc = min(Fs_in, Fs_out) / 2; % Fréquence de coupure
N = 128; % Ordre du filtre
lpFilt = designfilt('lowpassfir', 'PassbandFrequency', Fc - 2000, ...
    'StopbandFrequency', Fc + 2000, 'SampleRate', Fs_in * L, 'FilterOrder', N);

% Appliquer le filtre passe-bas
signal_filtered = filter(lpFilt, signal_upsampled);

% Sous-échantillonnage (garder 1 échantillon sur M)
signal_downsampled = downsample(signal_filtered, M);

%% 2. Version améliorée
% Conception d'un SRC avec approche polyphasique pour optimiser les calculs
P = L; % Facteur d'interpolation
Q = M; % Facteur de décimation

% Appliquer un SRC polyphasique
signal_polyphase = resample(input_signal, P, Q);

% Comparaison des résultats des deux méthodes
N_comp = min(length(signal_downsampled), length(signal_polyphase)); % Taille de comparaison
comparison_naive_polyphase = norm(signal_downsampled(1:N_comp) - signal_polyphase(1:N_comp)) / N_comp;

%% 3. Version optimale
% Pour l'implémentation optimisée sur cible, remplacer le filtre FIR par une implémentation directe
% Ajuster l'algorithme en utilisant des étapes prêtes pour une intégration matérielle
% Par exemple, utiliser `fir1` avec des coefficients en nombres fixes
N_opt = 64; % Ordre réduit pour minimiser la complexité
h_opt = fir1(N_opt, Fc / (Fs_in)); % FIR optimal pour intégration simple

% Sur-échantillonnage par interpolation (zero-stuffing)
interpolated = upsample(input_signal, P);

% Appliquer le filtre optimisé
filtered_opt = conv(interpolated, h_opt, 'same');

% Sous-échantillonnage (downsample par Q)
signal_optimized = filtered_opt(1:Q:end);

%% Comparaison des résultats
figure;
subplot(2,1,1);
plot(signal_polyphase, 'b'); hold on;
plot(signal_optimized, '--r');
xlabel('Samples');
ylabel('Amplitude');
legend('Polyphasique', 'Optimisée');
title('Comparaison des sorties SRC');

% Affichage des spectres
subplot(2,1,2);
fft_poly = fftshift(abs(fft(signal_polyphase))); 
fft_opt = fftshift(abs(fft(signal_optimized)));
plot(linspace(-Fs_out/2, Fs_out/2, length(fft_poly)), 20*log10(fft_poly), 'b'); hold on;
plot(linspace(-Fs_out/2, Fs_out/2, length(fft_opt)), 20*log10(fft_opt), '--r');
xlabel('Fréquence (Hz)');
ylabel('Amplitude (dB)');
legend('Polyphasique', 'Optimisée');
title('Spectres des sorties SRC');

%% Évaluation qualité des méthodes
% Calcul SNR pour chaque méthode
snr_naive = snr(signal_downsampled, data.w441(50e3:50e3+N_comp));
snr_polyphase = snr(signal_polyphase, data.w441(50e3:50e3+N_comp));
snr_opt = snr(signal_optimized, data.w441(50e3:50e3+N_comp));

disp('--- Évaluation des solutions ---');
fprintf('Naïve SNR : %.2f dB\n', snr_naive);
fprintf('Polyphasique SNR : %.2f dB\n', snr_polyphase);
fprintf('Optimale SNR : %.2f dB\n', snr_opt);
