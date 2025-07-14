%% itd_and_ild_bode_plots.m


small_itd_delay = 50e-6; % seconds
large_itd_delay = 500e-6; % seconds

Fs = 44100;
%
click_input = zeros(1, Fs);
click_output = zeros(1, Fs);
click_input(1) = 1;
delay_samples = round(small_itd_delay * Fs);
click_output(1 + delay_samples) = 1;
N = length(click_input);
nfft = 2^nextpow2(N); %  Find next power of 2 for FFT length
fft_input = fft(click_input,nfft);
fft_output = fft(click_output,nfft);
f = Fs*(0:(nfft/2))/nfft;


magnitude_dB_input = 20*log10(abs(fft_input(1:nfft/2+1)));
phase_deg_input = unwrap(angle(fft_input(1:nfft/2+1)))*180/pi;

magnitude_dB_output = 20*log10(abs(fft_output(1:nfft/2+1)));
phase_deg_output = unwrap(angle(fft_output(1:nfft/2+1)))*180/pi;

figure;
% plot small ITDs
subplot(2,4,1)
semilogx(f,magnitude_dB_input - magnitude_dB_output,'LineWidth',2)
xlim([0,8000])
ylabel('Magnitude (dB)','FontSize',18)
title('Small ITDs','FontSize',18)
ylim([-5,40])

subplot(2,4,5)
semilogx(f,phase_deg_output - phase_deg_input,'LineWidth',2)
xlim([0,8000])
ylabel('Phase (deg)','FontSize',18)
ylim([-1000,100])



click_input = zeros(1, Fs);
click_output = zeros(1, Fs);
click_input(1) = 1;
delay_samples = round(large_itd_delay * Fs);
click_output(1 + delay_samples) = 1;
N = length(click_input);
nfft = 2^nextpow2(N); %  Find next power of 2 for FFT length
fft_input = fft(click_input,nfft);
fft_output = fft(click_output,nfft);
f = Fs*(0:(nfft/2))/nfft;


magnitude_dB_input = 20*log10(abs(fft_input(1:nfft/2+1)));
phase_deg_input = unwrap(angle(fft_input(1:nfft/2+1)))*180/pi;

magnitude_dB_output = 20*log10(abs(fft_output(1:nfft/2+1)));
phase_deg_output = unwrap(angle(fft_output(1:nfft/2+1)))*180/pi;
% plot large ITDs
subplot(2,4,2)
semilogx(f,magnitude_dB_input - magnitude_dB_output,'LineWidth',2)
xlim([0,8000])
title('Large ITDs','FontSize',18)
ylim([-5,40])

subplot(2,4,6)
semilogx(f,phase_deg_output - phase_deg_input,'LineWidth',2)
xlim([0,8000])
ylim([-1000,100])


% Natural ILDs
hrtf_data = table2array(readtable('C:\Users\benri\Documents\GitHub\MILD-Master\stim\hrtf_kemar_0el.csv'));
hrtf_r_ind = mod(mod(70,360)/5,(360/5));
hrtf_l_ind = mod(360 - mod(70,360)/5,360/5);

% this HRTF data is the impulse response of the syhstem (I.e., we convolve
% it)
hrtf_l_data = hrtf_data(hrtf_l_ind,:);
hrtf_r_data = hrtf_data(hrtf_r_ind,:);

% make a click
impulse_response_L = hrtf_l_data;
impulse_response_R = hrtf_r_data;

N = length(impulse_response_L);
nfft = 2^nextpow2(N); %  Find next power of 2 for FFT length
H_f_L = fft(impulse_response_L, nfft);
H_f_R = fft(impulse_response_R, nfft);

f = Fs*(0:(nfft/2))/nfft;
magnitude_dB_L = 20*log10(abs(H_f_L(1:nfft/2+1)));
phase_deg_L = unwrap(angle(H_f_L(1:nfft/2+1)))*180/pi;

magnitude_dB_R = 20*log10(abs(H_f_R(1:nfft/2+1)));
phase_deg_R= unwrap(angle(H_f_R(1:nfft/2+1)))*180/pi;

sys_frd_L = frd(H_f_L(1:nfft/2+1), 2*pi*f); 
sys_frd_R = frd(H_f_R(1:nfft/2+1), 2*pi*f); 


subplot(2,4,3)
semilogx(f,magnitude_dB_L - magnitude_dB_R,'LineWidth',2)
xlim([0,8000])
title('Natural ILDs','FontSize',18)

subplot(2,4,7)
semilogx(f,zeros(1,length(phase_deg_L)),'LineWidth',2) % phase_deg_L - phase_deg_R
xlim([0,8000])
ylim([-1000,100])


subplot(2,4,4)
semilogx(f,10*ones(1,length(magnitude_dB_L)),'LineWidth',2)
xlim([0,8000])
ylim([-5,40])
title('Broadband ILDs','FontSize',18)

subplot(2,4,8)
semilogx(f,zeros(1,length(phase_deg_L)),'LineWidth',2)
xlim([0,8000])
ylim([-1000,100])

