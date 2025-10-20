
[blue,~] = audioread('C:\Users\benri\Documents\GitHub\Broadband-ILD-fNIRS\unprocessed stim\unprocessed stim\bob_all\blue_short.wav');
[shoe,fs] = audioread('C:\Users\benri\Documents\GitHub\Broadband-ILD-fNIRS\unprocessed stim\unprocessed stim\bob_all\shoe_short.wav');

blue_lead = [blue; zeros(1,fs*0.15)'];
blue_lag = [zeros(1,fs*0.15)'; blue];
shoe_lead = [shoe; zeros(1,fs*0.15)'];
shoe_lag = [zeros(1,fs*0.15)'; shoe];
mask_thresh = -10;

win = round(fs/70);
noverlap = round(fs/140);
hop = win - noverlap;
L = length(blue);

nframes = floor((L - win)/hop) + 1;             % full frames only
L_trim = (nframes-1)*hop + win;                 % length that gives exact frames
x_trim = blue(1:L_trim);

[s_blue,f,t_blue] = spectrogram(x_trim,fs/70,round(fs/140),fs*10,fs,'onesided');
s_blue_mag = abs(s_blue);
s_blue_power = 10*log10(abs(s_blue_mag).^2);
s_blue_power(s_blue_power < mask_thresh) = min(s_blue_power,[],'all');

[s_shoe,f,t_shoe] = spectrogram(shoe,fs/70,round(fs/140),fs*10,fs,'onesided');
s_shoe_mag = abs(s_shoe);
s_shoe_power = 20*log10(s_shoe_mag);
s_shoe_power(s_shoe_power < mask_thresh) = min(s_shoe_power,[],'all');

% transform word blue power to RGB values
pmin = min(s_blue_power,[],'all');   % lowest expected power
pmax = max(s_blue_power,[],'all');    % highest expected power
blue_powerNorm = (s_blue_power - pmin) / (pmax - pmin);

% Just a little color
R = 1 - blue_powerNorm;              % fades from 1 → 0
G = 1 - blue_powerNorm;              % fades from 1 → 0
B = 0.2*ones(size(blue_powerNorm));      % always 1

blue_Image_little_color = cat(3, R, G, B);

R = 1 - blue_powerNorm;              % fades from 1 → 0
G = 1 - blue_powerNorm;              % fades from 1 → 0
B = ones(size(blue_powerNorm));      % always 1

blue_Image_full_color = cat(3, R, G, B);

natural_ilds = [linspace(0,1,find(f == 8000)),ones(1,length(f)-find(f==8000))]; % linearly increase for now up to 8 kHz
R = 1 - blue_powerNorm;              % fades from 1 → 0
G = 1 - blue_powerNorm;              % fades from 1 → 0
B = repmat(natural_ilds,size(blue_powerNorm,2),1)';

blue_Image_high_freq = cat(3, R, G, B);


%% transform word blue power to RGB values
pmin = min(s_blue_power,[],'all');
pmax = max(s_blue_power,[],'all');
blue_powerNorm = (s_blue_power - pmin) / (pmax - pmin);

% Target color for "blue" → (255,255,98)
target_R = 255/255; target_G = 255/255; target_B = 98/255;

% --- Just a little color (fades toward target color)
R = 0.3*blue_powerNorm*(target_R);
G = 0.3*blue_powerNorm*(target_G);
B = 0.3*blue_powerNorm*(target_B);
blue_Image_little_color = cat(3, R, G, B);

% --- Full vivid color (saturates to target color)
R = 1 - blue_powerNorm*(1 - target_R);
G = 1 - blue_powerNorm*(1 - target_G);
B = 1 - blue_powerNorm*(1 - target_B);
blue_Image_full_color = cat(3, R, G, B);

% --- Natural ILDs (higher freq coloration)
natural_ilds = [linspace(0,1,find(f == 8000)),ones(1,length(f)-find(f==8000))];
R = (1 - blue_powerNorm) + natural_ilds' * (target_R);
G = (1 - blue_powerNorm) + natural_ilds' * (target_G);
B = (1 - blue_powerNorm) + natural_ilds' * (target_B);
blue_Image_high_freq = cat(3, R, G, B);

%% transform word shoe power to RGB values
pmin = min(s_shoe_power,[],'all');
pmax = max(s_shoe_power,[],'all');
shoe_powerNorm = (s_shoe_power - pmin) / (pmax - pmin);

% Target color for "shoe" → (211,95,183)
target_R = 211/255; target_G = 95/255; target_B = 183/255;

% --- Just a little color (fades toward target color)
R = 0.4*shoe_powerNorm*(target_R);
G = 0.4*shoe_powerNorm*(target_G);
B = 0.4*shoe_powerNorm*(target_B);
shoe_Image_little_color = cat(3, R, G, B);

% --- Full vivid color
R = 1 - shoe_powerNorm*(1 - target_R);
G = 1 - shoe_powerNorm*(1 - target_G);
B = 1 - shoe_powerNorm*(1 - target_B);
shoe_Image_full_color = cat(3, R, G, B);

% --- Natural ILDs (higher freq coloration)
natural_ilds = [linspace(0,1,find(f == 5000)),ones(1,length(f)-find(f==5000))];
R = (1 - shoe_powerNorm) + natural_ilds' * (target_R);
G = (1 - shoe_powerNorm) + natural_ilds' * (target_G);
B = (1 - shoe_powerNorm) + natural_ilds' * (target_B);
shoe_Image_high_freq = cat(3, R, G, B);



% Plot each black on its own
fig = figure;
h1 = imagesc(t_blue,f,blue_powerNorm);
colormap(gray)
c = colormap(gray);
c_flipped = flipud(c);
colormap(c_flipped)
ax = gca;
ax.YDir = 'normal';
n = 256;
ylim([0,8000])
xlim([0,0.300])
set(h1, 'AlphaData', blue_powerNorm);
title('blue')

fig = figure;
h1 = imagesc(t_shoe,f,shoe_powerNorm);
colormap(gray)
c = colormap(gray);
c_flipped = flipud(c);
colormap(c_flipped)
ax = gca;
ax.YDir = 'normal';
n = 256;
ylim([0,8000])
xlim([0,0.300])
set(h1, 'AlphaData', shoe_powerNorm);
title('shoe')

% Plot both with JUST A LITTLE BIT COLOR 
fig = figure;
h1 = imagesc(t_blue,f,blue_Image_little_color);
colormap(gray)
hold on
h2 = imagesc(t_blue + 0.150,f,shoe_Image_little_color);
c = colormap(gray);
c_flipped = flipud(c);
colormap(c_flipped)
ax = gca;
ax.YDir = 'normal';
n = 256;
ylim([0,8000])
xlim([0,0.450])
set(h1, 'AlphaData', blue_powerNorm);
set(h2, 'AlphaData', shoe_powerNorm);


% Plot where color is greater with frequency (natural ILDs)
fig = figure;
h1 = imagesc(t_blue,f,blue_Image_high_freq);
hold on
h2 = imagesc(t_blue + 0.150,f,shoe_Image_high_freq);
ax = gca;
ax.YDir = 'normal';
n = 256;
ylim([0,8000])
xlim([0,0.450])
set(h1, 'AlphaData', blue_powerNorm);
set(h2, 'AlphaData', shoe_powerNorm);

% Plot one word blue and one red (Large ITDs/Broadband ILDs)
fig = figure;
h1 = imagesc(t_blue,f,blue_Image_full_color);
hold on
h2 = imagesc(t_blue + 0.150,f,shoe_Image_full_color);
ax = gca;
ax.YDir = 'normal';
n = 256;
ylim([0,8000])
xlim([0,0.450])
set(h1, 'AlphaData', blue_powerNorm);
set(h2, 'AlphaData', shoe_powerNorm);
% alpha('scaled')
% alphamap(linspace(0,1,n)); % transparency ramp
