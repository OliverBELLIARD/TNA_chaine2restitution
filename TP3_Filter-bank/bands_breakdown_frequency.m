clear all; close all; clc;

% Parameters
x = sin(2*pi*0.01*(0:999))'; % Test signal
M = 256; % Window size
% R = M/2; % 50% Overlap
R = M/4; % 75% Overlap
w_analysis = hamming(M); % Analysis window
w_synthesis = w_analysis.^2; % Example synthesis window: squared version

%% STFT Analysis
num_frames = floor((length(x) - M) / R) + 1;
X = zeros(M, num_frames);
for k = 0:num_frames-1
    x_seg = x(k*R + (1:M)); % Extract frame
    x_seg_windowed = x_seg .* w_analysis; % Apply analysis window
    X(:,k+1) = fft(x_seg_windowed); % FFT
end

% Signal processing here

%% ISTFT and WOLA Reconstruction
x_reconstructed = zeros(length(x), 1); % Initialize
for k = 0:num_frames-1
    x_ifft = ifft(X(:,k+1)); % Inverse FFT
    x_overlap = real(x_ifft) .* w_synthesis; % Apply synthesis window
    x_reconstructed(k*R + (1:M)) = x_reconstructed(k*R + (1:M)) + x_overlap;
end

%% Plot results
figure;
plot(x, 'b'); hold on;
plot(x_reconstructed, 'r--');
title('WOLA Reconstruction vs Original');
legend('Original Signal', 'Reconstructed Signal');
xlabel('Samples');
ylabel('Amplitude');