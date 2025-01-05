clear all; clc; close all;

%% Load Signal
file_path = 'pcm_48k.mat'; % Adjust to the correct file path
if exist(file_path, 'file') == 2
    load(file_path, '-mat'); % Load the file
else
    error('File not found: %s', file_path);
end

if ~exist('pcm_48', 'var')
    error('Variable pcm_48 does not exist in the loaded file.');
end

Xin = pcm_48(5e4:1e5); % Signal
Fs = 48000;   % Sampling frequency
N = length(Xin);
t = (0:N-1) / Fs;

%% Frequencies for octave bands (20 Hz to 20 kHz)
f_min = 20;                   % Minimum frequency
num_bands = 10;               % Number of bands
frequencies = f_min * 2.^(0:(num_bands-1)); % Octave band frequencies

%% Design Parameters
A = 60; % Stopband attenuation in dB
transition_width = 500; % Transition width in Hz
beta = 0.1102 * (A - 8.7); % Kaiser beta

% Dynamic order calculation for FIR
order = ceil((A - 7.95) / (2.285 * (transition_width / Fs))); % Dynamic order

% Filters
filters_lp = cell(1, num_bands); % Low-pass filters
filters_hp = [];                % High-pass filter for the last band

for i = 1:num_bands
    if i == num_bands
        % Low-pass filter for the last band
        cutoff_lp = frequencies(i) / (Fs / 2); % Normalized cutoff
        filters_lp{i} = fir1(order, cutoff_lp, 'low', kaiser(order + 1, beta)); % Low-pass filter

        % High-pass filter for the highest frequencies
        cutoff_hp = frequencies(i) / (Fs / 2); % Same cutoff for high-pass
        filters_hp = fir1(order, cutoff_hp, 'high', kaiser(order + 1, beta)); % High-pass filter
    else
        % Low-pass filter for intermediate bands
        cutoff = frequencies(i) / (Fs / 2); % Normalized cutoff
        filters_lp{i} = fir1(order, cutoff, 'low', kaiser(order + 1, beta)); % Low-pass filter
    end
end

%% Apply filters and calculate bands
subbands = cell(1, num_bands);
for i = 1:num_bands
    if i == 1
        % First band: Low-pass only
        filtered_upper = filter(filters_lp{i}, 1, Xin);
        subbands{i} = circshift(filtered_upper, ceil(-order / 2));
    elseif i == num_bands
        % Last band: Combine low-pass and high-pass
        filtered_lower = filter(filters_lp{i}, 1, Xin); % Low-pass for transition
        filtered_upper = filter(filters_hp, 1, Xin);    % High-pass for highest content
        subbands{i} = circshift(filtered_lower, ceil(-order / 2))...
        + circshift(filtered_upper, ceil(-order / 2));
    else
        % Intermediate bands: Subtract adjacent low-pass filtered signals
        filtered_upper = filter(filters_lp{i}, 1, Xin);
        filtered_lower = filter(filters_lp{i-1}, 1, Xin);
        subbands{i} = circshift(filtered_upper, ceil(-order / 2))...
        - circshift(filtered_lower, ceil(-order / 2));
    end
end

% Reconstruction from bands
y_reconstructed = zeros(size(Xin));
for i = 1:num_bands
    y_reconstructed = y_reconstructed + subbands{i};
end

%% Visualization of Results
figure;

% (1) Original vs Reconstructed Signal
subplot(3,1,1);
plot(t, Xin, 'b', 'LineWidth', 1.5); hold on;
plot(t, y_reconstructed, 'r', 'LineWidth', 1.5);
legend('Original Signal', 'Reconstructed Signal');
title('Original vs Reconstructed Signal');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 500/Fs]);
grid on;

% (2) Reconstruction Error
subplot(3,1,2);
error_signal = Xin - y_reconstructed;
plot(t, error_signal, 'k', 'LineWidth', 1.5);
title('Reconstruction Error');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% (3) Band Filters' Responses
subplot(3,1,3);
hold on;
for i = 1:num_bands-1
    [H, f] = freqz(filters_lp{i}, 1, 4096, Fs);
    semilogx(f, 20*log10(abs(H)), 'LineWidth', 1.5);
end
% Add responses for last band's filters
[H_hp, f] = freqz(filters_hp, 1, 4096, Fs);
semilogx(f, 20*log10(abs(H_hp)), 'LineWidth', 1.5);
[H_lp_last, f] = freqz(filters_lp{num_bands}, 1, 4096, Fs);
semilogx(f, 20*log10(abs(H_lp_last)), 'LineWidth', 1.5);

title('Bode Diagram of Kaiser FIR Filters');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([10, Fs/2]);
grid on;
legend([arrayfun(@(x) sprintf('Band %d (LP)', x), 1:num_bands-1, 'UniformOutput', false), ...
        {'Band 10 (HP)', 'Band 10 (LP)'}]);
hold off;

