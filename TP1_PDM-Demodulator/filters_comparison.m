clear all; clc; close all;

%% Load filters
fir1 = FIR_equiripple_5680();
fir2 = FIR_constrained_equiripple_5000();
fir3 = FIR_generalized_equiripple_11000();

iir1 = IIR_Butterworth_94();
iir2 = IIR_Chebyshev1_29();
iir3 = IIR_Chebyshev2_29();
iir4 = IIR_Eliptic_15();

filters_fir = {fir1, fir2, fir3}; % FIR filters
filters_iir = {iir1, iir2, iir3, iir4}; % IIR filters
labels_fir = {'FIR Equiripple 5680', 'FIR Constrained Equiripple 5000', 'FIR Generalized Equiripple 11000'};
labels_iir = {'IIR Butterworth 94', 'IIR Chebyshev1 29', 'IIR Chebyshev2 29', 'IIR Elliptic 15'};

Fs = 6.144e6; % Sampling frequency

%% Plot FIR filters: Gain and Phase
figure('Name', 'FIR Filters: Gain and Phase');
subplot(2, 1, 1);
hold on;
for i = 1:length(filters_fir)
    [H, f] = freqz(filters_fir{i}.Numerator, 1, 2048, Fs);
    plot(f / 1e3, 20 * log10(abs(H)), 'DisplayName', labels_fir{i});
end
xlabel('Frequency (kHz)');
ylabel('Gain (dB)');
title('FIR Filters: Gain');
grid on;
legend;
xlim([0 50]);
ylim([-150 10]);
hold off;

subplot(2, 1, 2);
hold on;
for i = 1:length(filters_fir)
    [H, f] = freqz(filters_fir{i}.Numerator, 1, 2048, Fs);
    plot(f / 1e3, unwrap(angle(H)) * (180 / pi), 'DisplayName', labels_fir{i});
end
xlabel('Frequency (kHz)');
ylabel('Phase (Degrees)');
title('FIR Filters: Phase');
grid on;
legend;
xlim([0 50]);
hold off;

%% Plot IIR filters: Gain and Phase
figure('Name', 'IIR Filters: Gain and Phase');
subplot(2, 1, 1);
hold on;
for i = 1:length(filters_iir)
    [H, f] = freqz(filters_iir{i});
    plot(f / 1e-3, 20 * log10(abs(H)), 'DisplayName', labels_iir{i});
end
xlabel('Frequency (kHz)');
ylabel('Gain (dB)');
title('IIR Filters: Gain');
grid on;
legend;
xlim([0 50]);
ylim([-150 10]);
hold off;

subplot(2, 1, 2);
hold on;
for i = 1:length(filters_iir)
    [H, f] = freqz(filters_iir{i});
    plot(f / 1e-3, unwrap(angle(H)) * (180 / pi), 'DisplayName', labels_iir{i});
end
xlabel('Frequency (kHz)');
ylabel('Phase (Degrees)');
title('IIR Filters: Phase');
grid on;
legend;
xlim([0 50]);
hold off;

%% Compare all filters together
figure('Name', 'All Filters: Gain and Phase');
subplot(2, 1, 1);
hold on;
for i = 1:length(filters_fir)
    [H, f] = freqz(filters_fir{i}.Numerator, 1, 2048, Fs);
    plot(f / 1e3, 20 * log10(abs(H)), 'DisplayName', labels_fir{i});
end
for i = 1:length(filters_iir)
    [H, f] = freqz(filters_iir{i});
    plot(f / 1e-3, 20 * log10(abs(H)), '--', 'DisplayName', labels_iir{i});
end
xlabel('Frequency (kHz)');
ylabel('Gain (dB)');
title('All Filters: Gain');
grid on;
legend;
xlim([0 50]);
ylim([-150 10]);
hold off;

subplot(2, 1, 2);
hold on;
for i = 1:length(filters_fir)
    [H, f] = freqz(filters_fir{i}.Numerator, 1, 2048, Fs);
    plot(f / 1e3, unwrap(angle(H)) * (180 / pi), 'DisplayName', labels_fir{i});
end
for i = 1:length(filters_iir)
    [H, f] = freqz(filters_iir{i});
    plot(f / 1e-3, unwrap(angle(H)) * (180 / pi), '--', 'DisplayName', labels_iir{i});
end
xlabel('Frequency (kHz)');
ylabel('Phase (Degrees)');
title('All Filters: Phase');
grid on;
legend;
xlim([0 50]);
hold off;
