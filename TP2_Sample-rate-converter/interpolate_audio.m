function Xout = interpolate_audio(Xin, Fs, L, Nbits)
%% Interpolate an audio signal by an L ratio with an SNR of Nbits resolution
if nargin < 4
    Nbits = 16; % Default: 16 bits for calculating stopband attenuation
end

% Upsampled sampling rate
Fsup = Fs * L;

%% Anti-Aliasing Filter Design
Apass = 0.01;        % Passband ripple in dB
Astop = 1.76 + 6.02 * Nbits; % Stopband attenuation in dB

Fpass = 20000;             % Passband Frequency
Fstop = Fs / 2;            % Stopband Frequency
Dpass = power(10, -Apass/20);  % Passband Ripple
Dstop = power(10, -Astop/20);  % Stopband Attenuation
dens  = 20;                 % Density Factor

% Calculate the filter order and coefficients
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fsup/2), [1 0], [Dpass, Dstop])
b = firpm(N, Fo, Ao, W, {dens}); % FIR filter coefficients

%% Upsample and Filter
% Zero-padding for upsampling
upsampled_signal = zeros(1, L * length(Xin));
upsampled_signal(1:L:end) = Xin; % Insert original samples with L-1 zeros in between

% Apply anti-aliasing filter
filtered_signal = filter(b, 1, upsampled_signal); % FIR filtering for anti-aliasing

%% Delay Compensation
group_delay = ceil(N / 2); % Delay introduced by the FIR filter
filtered_signal = filtered_signal(group_delay + 1:end); % Compensate for delay
filtered_signal = [filtered_signal, zeros(1, group_delay)]; % Align signal length

%% Gain Normalization to Preserve Amplitude
gain_compensation = max(filtered_signal) / max(Xin);  % Calculate the amplitude loss during filtering
Xout = filtered_signal / gain_compensation; % Apply gain compensation
end
