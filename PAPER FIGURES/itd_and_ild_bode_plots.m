%% itd_and_ild_bode_plots_fixed.m
% Fixed semilogx plotting for ITD and ILD analysis

clear; close all; clc;

%% Parameters
small_itd_delay = 50e-6;  % seconds
large_itd_delay = 500e-6; % seconds
Fs = 44100;               % sampling rate

%% --- SMALL ITD CASE -------------------------------------------------------
click_input = zeros(1, Fs);
click_output = zeros(1, Fs);
click_input(1) = 1;
delay_samples = round(small_itd_delay * Fs);
click_output(1 + delay_samples) = 1;

nfft = 2^14;
fft_input = fft(click_input, nfft);
fft_output = fft(click_output, nfft);
f = Fs * (0:(nfft/2)) / nfft;  % frequency vector (Hz)

magnitude_dB_input = 20*log10(abs(fft_input(1:nfft/2+1)));
phase_deg_input = unwrap(angle(fft_input(1:nfft/2+1))) * 180/pi;

magnitude_dB_output = 20*log10(abs(fft_output(1:nfft/2+1)));
phase_deg_output = unwrap(angle(fft_output(1:nfft/2+1))) * 180/pi;

% Skip DC bin (f=0) for plotting
f_plot = f(2:end);
mag_diff_small = magnitude_dB_input(2:end) - magnitude_dB_output(2:end);
phase_diff_small = phase_deg_output(2:end) - phase_deg_input(2:end);

figure;

% Magnitude
subplot(2,4,1)
semilogx(f_plot, mag_diff_small, 'LineWidth', 2, 'Color', 'k')
xlim([10 1000])
ylabel('Magnitude (dB)', 'FontSize', 18)
ylim([-5 30])
grid on

% Phase
subplot(2,4,5)
semilogx(f_plot, phase_diff_small, 'LineWidth', 2, 'Color', 'k')
xlim([10 1000])
ylabel('Phase (deg)', 'FontSize', 18)
ylim([-175 50])
grid on


%% --- LARGE ITD CASE -------------------------------------------------------
click_input = zeros(1, Fs);
click_output = zeros(1, Fs);
click_input(1) = 1;
delay_samples = round(large_itd_delay * Fs);
click_output(1 + delay_samples) = 1;

fft_input = fft(click_input, nfft);
fft_output = fft(click_output, nfft);

magnitude_dB_input = 20*log10(abs(fft_input(1:nfft/2+1)));
phase_deg_input = unwrap(angle(fft_input(1:nfft/2+1))) * 180/pi;
magnitude_dB_output = 20*log10(abs(fft_output(1:nfft/2+1)));
phase_deg_output = unwrap(angle(fft_output(1:nfft/2+1))) * 180/pi;

mag_diff_large = magnitude_dB_input(2:end) - magnitude_dB_output(2:end);
phase_diff_large = phase_deg_output(2:end) - phase_deg_input(2:end);

% Magnitude
subplot(2,4,2)
semilogx(f_plot, mag_diff_large, 'LineWidth', 2, 'Color', 'k')
xlim([10 1000])
ylim([-5 30])
grid on

% Phase
subplot(2,4,6)
semilogx(f_plot, phase_diff_large, 'LineWidth', 2, 'Color', 'k')
xlim([10 1000])
ylim([-175 50])
grid on


%% --- NATURAL ILDs ---------------------------------------------------------
hrtf_data = table2array(readtable('C:\Users\benri\Documents\GitHub\MILD-Master\stim\hrtf_kemar_0el.csv'));
hrtf_r_ind = mod(mod(70,360)/5, (360/5));
hrtf_l_ind = mod(360 - mod(70,360)/5, 360/5);

hrtf_l_data = hrtf_data(hrtf_l_ind,:);
hrtf_r_data = hrtf_data(hrtf_r_ind,:);

impulse_response_L = hrtf_l_data;
impulse_response_R = hrtf_r_data;

N = length(impulse_response_L);
nfft = 2^14; % 2^nextpow2(N);
H_f_L = fft(impulse_response_L, nfft);
H_f_R = fft(impulse_response_R, nfft);

f = Fs * (0:(nfft/2)) / nfft;
f_plot = f(2:end);

magnitude_dB_L = 20*log10(abs(H_f_L(1:nfft/2+1)));
phase_deg_L = zeros(1,nfft/2+1); % unwrap(angle(H_f_L(1:nfft/2+1))) * 180/pi;
magnitude_dB_R = 20*log10(abs(H_f_R(1:nfft/2+1)));
phase_deg_R = zeros(1,nfft/2+1); % unwrap(angle(H_f_R(1:nfft/2+1))) * 180/pi;

mag_diff_ild = magnitude_dB_L(2:end) - magnitude_dB_R(2:end);
phase_diff_ild = phase_deg_L(2:end) - phase_deg_R(2:end);

% Magnitude
subplot(2,4,3)
semilogx(f_plot, mag_diff_ild, 'LineWidth', 2, 'Color', 'k')
xlim([10 8000])
ylim([-5 30])
grid on

% Phase
subplot(2,4,7)
semilogx(f_plot, phase_diff_ild, 'LineWidth', 2, 'Color', 'k')
xlim([10 8000])
ylim([-175 50])
grid on


%% --- BROADBAND ILDs -------------------------------------------------------
broadband_mag = 10 * ones(size(f_plot));
broadband_phase = zeros(size(f_plot));

% Magnitude
subplot(2,4,4)
semilogx(f_plot, broadband_mag, 'LineWidth', 2, 'Color', 'k')
xlim([10 8000])
ylim([-5 30])
grid on

% Phase
subplot(2,4,8)
semilogx(f_plot, broadband_phase, 'LineWidth', 2, 'Color', 'k')
xlim([10 8000])
ylim([-175 50])
grid on
