%% Load audio
[blue,~] = audioread('/Users/benri/Documents/GitHub/Broadband-ILD-fNIRS/unprocessed stim/unprocessed stim/bob_all/shoe_short.wav');
[shoe,fs] = audioread('/Users/benri/Documents/GitHub/Broadband-ILD-fNIRS/unprocessed stim/unprocessed stim/bob_all/blue_short.wav');

mask_thresh = -15;

%% Spectrogram parameters
win = round(fs/70);
noverlap = round(fs/140);
hop = win - noverlap;
L = length(blue);
nframes = floor((L - win)/hop) + 1;
L_trim = (nframes-1)*hop + win;
x_trim = blue(1:L_trim);

%% BLUE spectrogram
[s_blue,f,t_blue] = spectrogram(x_trim,fs/70,round(fs/140),fs*10,fs,'onesided');
s_blue_power = 10*log10(abs(s_blue).^2);
s_blue_power(s_blue_power < mask_thresh) = mask_thresh;

%% SHOE spectrogram
[s_shoe,f,t_shoe] = spectrogram(shoe,fs/70,round(fs/140),fs*10,fs,'onesided');
s_shoe_power = 20*log10(abs(s_shoe));
s_shoe_power(s_shoe_power < mask_thresh) = mask_thresh;

%% Normalize
blue_powerNorm = normalizeMatrix(s_blue_power);
shoe_powerNorm = normalizeMatrix(s_shoe_power);

%% Target colors
color_blue = [1.0, 1.0, 98/255];      % yellow
color_shoe = [211/255, 95/255, 183/255];  % magenta

%% Make color images
blue_Image_little_color = makeColorImageScaled(blue_powerNorm, color_blue, 0.5);
blue_Image_full_color   = makeColorImageScaled(blue_powerNorm, color_blue, 1.0);
blue_Image_high_freq    = makeILDImage(blue_powerNorm, color_blue, f, 8000);

shoe_Image_little_color = makeColorImageScaled(shoe_powerNorm, color_shoe, 0.5);
shoe_Image_full_color   = makeColorImageScaled(shoe_powerNorm, color_shoe, 1.0);
shoe_Image_high_freq    = makeILDImage(shoe_powerNorm, color_shoe, f, 5000);

%% Overlap mask
t_shift = 0.150; % sec
dt = t_blue(2)-t_blue(1);
shift_bins = round(t_shift / dt);
shoe_shifted = circshift(shoe_powerNorm, [0 shift_bins]);
shoe_shifted(:,1:shift_bins) = 0;
overlap_mask = (blue_powerNorm > 0.3) & (shoe_shifted > 0.3);
gray_val = 0.5;

%% Low-power background mask (both spectrograms near zero)
bg_mask = (blue_powerNorm < 0.05) & (shoe_shifted < 0.05);

%% --- Plot 1: LITTLE COLOR ---
img_little = overlay_gray_white(blue_Image_little_color, circshift(shoe_Image_little_color,[0 shift_bins]), ...
    overlap_mask, gray_val, bg_mask);
figure;
imagesc(t_blue,f,img_little);
set(gca,'YDir','normal');
ylim([0 8000]); xlim([0 0.45]);
title('Little Color (Gray Overlap, White Background)');
xlabel('Time (s)'); ylabel('Frequency (Hz)');

%% --- Plot 2: HIGH-FREQUENCY ILDs ---
img_high = overlay_gray_white(blue_Image_high_freq, circshift(shoe_Image_high_freq,[0 shift_bins]), ...
    overlap_mask, gray_val, bg_mask);
figure;
imagesc(t_blue,f,img_high);
set(gca,'YDir','normal');
ylim([0 8000]); xlim([0 0.45]);
title('High-Frequency ILDs (Gray Overlap, White Background)');
xlabel('Time (s)'); ylabel('Frequency (Hz)');

%% --- Plot 3: FULL VIVID COLOR ---
img_full = overlay_gray_white(blue_Image_full_color, circshift(shoe_Image_full_color,[0 shift_bins]), ...
    overlap_mask, gray_val, bg_mask);
figure;
imagesc(t_blue,f,img_full);
set(gca,'YDir','normal');
ylim([0 8000]); xlim([0 0.45]);
title('Full Vivid Color (Gray Overlap, White Background)');
xlabel('Time (s)'); ylabel('Frequency (Hz)');

%% ===== Helper Functions =====
function Mnorm = normalizeMatrix(M)
    Mnorm = (M - min(M(:))) / (max(M(:)) - min(M(:)));
end

function img = makeColorImageScaled(powerNorm, baseColor, intensityScale)
    R = baseColor(1) * powerNorm * intensityScale;
    G = baseColor(2) * powerNorm * intensityScale;
    B = baseColor(3) * powerNorm * intensityScale;
    img = cat(3,R,G,B);
end

function img = makeILDImage(powerNorm, baseColor, f, f_cut)
    idx = find(f <= f_cut,1,'last');
    if isempty(idx), idx = length(f); end
    natural_ilds = [linspace(0,1,idx), ones(1,length(f)-idx)];
    R = baseColor(1)*powerNorm .* natural_ilds';
    G = baseColor(2)*powerNorm .* natural_ilds';
    B = baseColor(3)*powerNorm .* natural_ilds';
    img = cat(3,R,G,B);
end

function img_out = overlay_gray_white(base1, base2, mask, gray_val, bg_mask)
    % Combine two RGB images
    img_out = max(base1, base2);
    % Set overlap pixels to gray
    for c = 1:3
        ch = img_out(:,:,c);
        ch(mask) = gray_val;
        img_out(:,:,c) = ch;
    end
    % Set low-power background to white
    for c = 1:3
        ch = img_out(:,:,c);
        ch(bg_mask) = 1.0;
        img_out(:,:,c) = ch;
    end
end
