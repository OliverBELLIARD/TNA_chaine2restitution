clear all; close all; clc;

%% Parameters
fs = 48000; % Sampling frequency
f_min = 20; % Minimum frequency for octave bands
num_bands = 10; % Number of octave bands
frequencies = f_min * 2.^(0:num_bands-1); % Center frequencies for each band

x = load('pcm_48k.mat', '-mat').pcm_48(5e4:1e5); % Signal input
M = 256; % Window size
overlap_fraction = 0.67; % 60% overlap
R = floor(M * (1 - overlap_fraction)); % Hop size

%% Window Functions
w_analysis = hamming(M); % Analysis window
w_synthesis = w_analysis.^2; % Synthesis window (squared version of analysis window)

%% STFT Analysis
num_frames = floor((length(x) - M) / R) + 1;
X = zeros(M, num_frames);
for k = 0:num_frames-1
    x_seg = x(k*R + (1:M)); % Extract frame
    x_seg_windowed = x_seg .* w_analysis; % Apply analysis window
    X(:,k+1) = fft(x_seg_windowed); % FFT
end

%% Frequency decomposition into 10 octave bands
H = zeros(M, num_bands); % Filters in the frequency domain
f = (0:M-1) * (fs / M); % Frequency axis
for i = 1:num_bands
    if i == 1
        % Low-pass filter for first band
        H(:,i) = (f <= frequencies(i+1));
    elseif i == num_bands
        % High-pass filter for last band
        H(:,i) = (f >= frequencies(i));
    else
        % Band-pass filter for intermediate bands
        H(:,i) = (f >= frequencies(i) & f < frequencies(i+1));
    end
end

%% Process each band
X_bands = cell(1, num_bands); % Decomposed bands in STFT domain
for i = 1:num_bands
    X_bands{i} = X .* H(:,i); % Apply frequency-domain filters
end

%% Reconstruct each band
x_bands = cell(1, num_bands);
for i = 1:num_bands
    x_reconstructed = zeros(length(x), 1); % Initialize for this band
    for k = 0:num_frames-1
        x_ifft = ifft(X_bands{i}(:,k+1)); % Inverse FFT
        x_overlap = real(x_ifft) .* w_synthesis; % Apply synthesis window
        x_reconstructed(k*R + (1:M)) = x_reconstructed(k*R + (1:M)) + x_overlap;
    end
    x_bands{i} = x_reconstructed; % Store the reconstructed signal for this band
end

%% Combine all bands for final reconstruction
x_final = zeros(size(x));
for i = 1:num_bands
    x_final = x_final + x_bands{i}; % Sum all bands
end

%% Plot Results

% (1) Plot each band in a single figure
figure;
num_rows = ceil(num_bands / 2); % Calculate the number of rows needed for 2 columns
for i = 1:num_bands
    subplot(num_rows, 2, i); % Position subplot in a grid of num_rows x 2
    plot(x_bands{i}); % Plot the sub-band
    title(sprintf('Band %d: %.1f Hz - %.1f Hz', i, frequencies(max(1,i-1)), frequencies(min(end,i+1))));
    xlabel('Samples');
    ylabel('Amplitude');
    xlim([0 50000])
    grid on;
end

% (2) Plot original signal, reconstructed signal, and reconstruction error
figure;
subplot(2,1,1);
plot(x, 'b'); hold on;
plot(x_final, 'r');
title('Original vs Reconstructed Signal');
legend('Original Signal', 'Reconstructed Signal');
xlabel('Samples');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(x - x_final, 'k'); % Error
title('Reconstruction Error');
xlabel('Samples');
ylabel('Amplitude');
grid on;

% (3) Plot frequency spectra of original and reconstructed signal
N_fft = 4096;
X_spectrum = abs(fft(x, N_fft)); % Spectrum of original
X_final_spectrum = abs(fft(x_final, N_fft)); % Spectrum of reconstructed

f_spectrum = linspace(0, fs/2, N_fft/2); % Frequency axis for the first half of FFT

figure;
semilogx(f_spectrum, 20*log10(X_spectrum(1:N_fft/2)), 'b', 'LineWidth', 1.5); hold on;
semilogx(f_spectrum, 20*log10(X_final_spectrum(1:N_fft/2)), 'r', 'LineWidth', 1.5);
title('Spectrum of Original vs Reconstructed Signal');
legend('Original Spectrum', 'Reconstructed Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
xlim([20, fs/2]);


