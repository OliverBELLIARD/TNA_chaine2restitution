clear all; clc; close all;

%% Load Signal
Xin = load("pcm_48k.mat", '-mat').pcm_48;% Input PCM signal
Fs = 48000;                              % Sampling frequency
N = length(Xin);                         % Signal length
t = (0:N-1)/Fs;                          % Time vector

%% Define Octave Bands (in Hz)
% Define cut-off frequencies for 6 bands
cutoff_freqs = [0, 375, 750, 1500, 3000, 6000, 12000]; % Band edges
num_bands = length(cutoff_freqs) - 1;

%% Design Filters
% Filter parameters
filter_order = 128; % Higher order for sharper transitions

% Create filters for each band
filters = cell(1, num_bands);
for k = 1:num_bands
    % Bandpass filter design
    low = cutoff_freqs(k) / (Fs/2); % Normalize by Nyquist frequency
    high = cutoff_freqs(k+1) / (Fs/2);
    if low == 0
        % Lowpass filter for the first band
        filters{k} = fir1(filter_order, high, 'low');
    elseif high == 1
        % Highpass filter for the last band
        filters{k} = fir1(filter_order, low, 'high');
    else
        % Bandpass filter for intermediate bands
        filters{k} = fir1(filter_order, [low, high]);
    end
end

%% Filter Signal into Bands
bands = cell(1, num_bands);
for k = 1:num_bands
    bands{k} = filter(filters{k}, 1, Xin); % Filter the signal
end

%% Reconstruct Signal from Bands
X_reconstructed = zeros(size(Xin));
for k = 1:num_bands
    X_reconstructed = X_reconstructed + bands{k};
end

%% Frequency Response of Filters
figure;
hold on;
N_fft = 2048;
f = linspace(0, Fs/2, N_fft/2); % Frequency axis
for k = 1:num_bands
    H = abs(fft(filters{k}, N_fft));
    semilogx(f, 20*log10(H(1:N_fft/2)), 'DisplayName', ['Band ', num2str(k)]);
end
title('Frequency Response of Octave Band Filters');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('show');
grid on;

%% Compare Original and Reconstructed Signals
figure;

subplot(3,1,1);
plot(t, Xin, 'b');
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(t, X_reconstructed, 'r');
title('Reconstructed Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,3);
plot(t, Xin - X_reconstructed, 'k');
title('Reconstruction Error');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% Display Reconstruction Error
error = max(abs(Xin - X_reconstructed));
disp(['Maximum Reconstruction Error: ', num2str(error)]);

%% --- Octave ---
function decompose_et_reconstruire_signal(nom_fichier)
  % Charger le signal audio
  [signal, fe] = audioread(nom_fichier);

  % Vérifier la fréquence d'échantillonnage
  if fe != 48000
    fprintf("La fréquence d'échantillonnage doit être de 48 kHz\n");
    return;
  endif

  % Créer les limites de fréquences des bandes d'octave
  f_min = 20; % Fréquence minimale audible
  f_max = fe / 2; % Fréquence maximale (limite de Nyquist)
  bandes = [];
  for i = 1:6
    bandes(i,:) = [f_min, f_min * 2];
    f_min = f_min * 2;
  endfor

  % Créer les filtres
  filtres = [];
  for i = 1:6
    filtres(i) = fir1(100, bandes(i,:) / (fe / 2), 'bandpass');
  endfor

  % Appliquer les filtres au signal
  signaux_filtres = [];
  for i = 1:6
    signaux_filtres(i,:) = filter(filtres(i), 1, signal);
  endfor

  % Reconstruire le signal
  signal_reconstruit = sum(signaux_filtres, 1);

  % Calculer l'erreur de reconstruction (erreur quadratique moyenne)
  erreur = mean((signal - signal_reconstruit).^2);

  % Visualiser les signaux
  figure;
  plot(signal);
  hold on;
  plot(signal_reconstruit, 'r');
  legend('Signal original', 'Signal reconstruit');
  title('Comparaison du signal original et reconstruit');
  xlabel('Échantillons');
  ylabel('Amplitude');

  fprintf("Erreur quadratique moyenne : %f\n", erreur);
endfunction

% Exemple d'utilisation
decompose_et_reconstruire_signal("mon_signal.wav");
