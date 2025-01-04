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

clear all; clc; close all;

% Paramètres généraux
Fs = 48000;  % Fréquence d'échantillonnage
N = 2^14;    % Taille du signal (puissance de 2 pour simplifier)
t = (0:N-1)/Fs;
Xin = cos(2*pi*1000*t) + cos(2*pi*5000*t) + cos(2*pi*15000*t); % Signal d'entrée (test)

% Conception des filtres
n = 50; % Ordre des filtres
lp = fir1(n, 0.5, 'low'); % Filtre passe-bas
hp = fir1(n, 0.5, 'high'); % Filtre passe-haut

% Fonction récursive pour appliquer la structure en arbre
function subbands = filter_tree(signal, lp, hp, levels)
    if levels == 0
        % Condition de terminaison : pas de découpe supplémentaire
        subbands = {signal};
    else
        % Filtrage passe-bas et passe-haut
        low_band = filter(lp, 1, signal);
        high_band = filter(hp, 1, signal);

        % Décimation par 2
        low_band = downsample(low_band, 2);
        high_band = downsample(high_band, 2);

        % Appel récursif pour les sous-bandes
        subbands_low = filter_tree(low_band, lp, hp, levels-1);
        subbands_high = filter_tree(high_band, lp, hp, levels-1);

        % Combiner les sous-bandes
        subbands = [subbands_low, subbands_high];
    end
end

% Nombre de niveaux d'arborescence
levels = 4;

% Application de la structure en arbre
subbands = filter_tree(Xin, lp, hp, levels);

% Reconstruction du signal
function y_reconstructed = reconstruct_tree(subbands, lp, hp)
    if length(subbands) == 1
        % Condition de terminaison : une seule bande
        y_reconstructed = subbands{1};
    else
        % Récupération des bandes passe-bas et passe-haut
        low_band = subbands(1:(length(subbands)/2));
        high_band = subbands((length(subbands)/2)+1:end);

        % Reconstruction des signaux
        y_low = reconstruct_tree(low_band, lp, hp);
        y_high = reconstruct_tree(high_band, lp, hp);

        % Suréchantillonnage par 2
        y_low_up = upsample(y_low, 2);
        y_high_up = upsample(y_high, 2);

        % Filtrage inverse
        y_low_filtered = filter(lp, 1, y_low_up);
        y_high_filtered = filter(hp, 1, y_high_up);

        % Combinaison
        y_reconstructed = y_low_filtered + y_high_filtered;
    end
end

% Reconstruction du signal
y_reconstructed = reconstruct_tree(subbands, lp, hp);

% Visualisation
figure;
subplot(3,1,1);
plot(t, Xin);
title('Signal d''entrée');
xlabel('Temps (s)'); ylabel('Amplitude');

subplot(3,1,2);
hold on;
for i = 1:length(subbands)
    plot(downsample(t, 2^(levels-1)), subbands{i}(1:downsample(end, 2^(levels-1))));
end
title('Signaux filtrés en sous-bandes');
xlabel('Temps (s)'); ylabel('Amplitude');

subplot(3,1,3);
plot(t, y_reconstructed);
title('Signal reconstruit');
xlabel('Temps (s)'); ylabel('Amplitude');

