%% Load Input Signal
clear all; clc; close all;

Xin = load("pcm_48k.mat", '-mat').pcm_48; % Input signal
Fs = 48000; % Sampling Frequency
N = length(Xin);

%% Define Octave Bands
% One octave bands: Center Frequencies (logarithmically spaced)
f_low = 20;        % Start frequency (Hz)
f_high = Fs/2;     % Nyquist frequency

% Generate center frequencies for octave bands
fc = []; 
while f_low < f_high
    fc = [fc f_low];
    f_low = f_low * 2; % Next octave
end

%% Filter Design
filter_order = 128; % FIR Filter order (choose even for linear phase)
Xbands = zeros(N, length(fc)); % Preallocate filtered signal matrix

% Iterate through octave bands and design FIR bandpass filters
for k = 1:length(fc)
    if k == 1 % Special case for the first band (lowpass)
        fpass = fc(k)*sqrt(2)/Fs; % Upper passband frequency (normalized)
        b = fir1(filter_order, fpass, 'low'); % Lowpass FIR filter
    elseif k == length(fc) % Last band (highpass)
        fpass = fc(k)/sqrt(2)/Fs; % Lower passband frequency (normalized)
        b = fir1(filter_order, fpass, 'high'); % Highpass FIR filter
    else % General case: Bandpass filter
        f1 = fc(k)/sqrt(2); % Lower edge of the band
        f2 = fc(k)*sqrt(2); % Upper edge of the band
        b = fir1(filter_order, [f1 f2]/(Fs/2), 'bandpass'); % Bandpass FIR
    end
    
    % Filter input signal to isolate band
    Xbands(:,k) = filter(b, 1, Xin); 
    
    % Optional: Display filter response
    [H, W] = freqz(b, 1, 2048, Fs);
    figure(1); hold on;
    loglog(W, abs(H), 'DisplayName', sprintf('Band %d: %.1f Hz', k, fc(k)));
end

%% Plot Filter Responses
title('Octave Filter Bank Frequency Responses');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on; legend show;

%% Visualize Results
% Time-domain plot for filtered bands
figure(2);
for k = 1:length(fc)
    subplot(length(fc),1,k);
    plot(Xbands(:,k));
    title(sprintf('Filtered Signal: Band %d (Center: %.1f Hz)', k, fc(k)));
    ylabel('Amplitude');
    grid on;
end
xlabel('Sample Index');

%% Perfect Reconstruction
% Sum the filtered outputs
Xreconstructed = sum(Xbands, 2);

% Verify reconstruction
figure(3);
subplot(2,1,1);
plot(Xin);
title('Original Signal');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(Xreconstructed);
title('Reconstructed Signal (Perfect Reconstruction)');
ylabel('Amplitude');
grid on;
