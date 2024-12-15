function Xout = decimate_audio(Xin, Fs, M, bits)
%% Filters and decimates an audio signal with a bit proportional SNR
if nargin < 4
    bits = 20;  % Defalut bits resolution for PCM
end

% Filter design
Fpass = 20000;    % Passband Frequency
Fstop = Fs/(2*M);    % Stopband Frequency
Apass = 0.01;     % Passband Ripple (dB)
Astop = 1.76 + 6.02*bits;   % Stopband Attenuation (dB), bit proportional SNR

h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop, Fs);

Hd = design(h, 'ellip', ...
    'MatchExactly', 'both', ...
    'SystemObject', true);

% lowpass Filtering
Xfilt = step(Hd, Xin);

% Decimation
Xout = decimate(Xfilt, M);