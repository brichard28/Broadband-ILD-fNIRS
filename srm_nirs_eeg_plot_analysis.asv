%% srm_nirs_eeg_plot_analysis

% Author: Benjamin Richardson
% Created 7/27/2023

% This script takes the epochs structure and group results structure from
% srm_nirs_eeg_analysis.py and plots results for the SRM-NIRS-EEG-1
% experiment

addpath C:\Users\benri\Documents\GitHub\Spatial-Hearing-fNIRS\errorbar_files

subject_ID = ["NDARCD778KPR","NDARLU426TBZ"]; %"NDARYZ656HJ9",
curr_subject_ID = ["NDARCD778KPR","NDARLU426TBZ"]; %"NDARYZ656HJ9",

%% Conditions Key
% '15.0':'Masker Noise ITD 50 Targ Left'
% '14.0':'Masker Noise ILD 10 Targ Right'
% '13.0':'Masker Noise ILD 10 Targ Left'
% '12.0':'Masker Speech ILD 10 Targ Right'
% '11.0':'Masker Speech ITD 500 Targ Left'
% '10.0':'Masker Noise ITD 500 Targ Left'
% '9.0':'Masker Speech ITD 50 Targ Left'
% '8.0':'Masker Noise ITD 500 Targ Right'
% '7.0':'Masker Speech ITD 50 Targ Right'
% '6.0':'Masker Speech ILD 10 Targ Left'
% '5.0':'Masker Speech ITD 500 Targ Right'
% '4.0':'Masker Noise ITD 50 Targ Right'
% '3.0':'Hold'
% '2.0':'Exhale'
% '1.0':'Inhale'
conditions = ["Inhale",...
"Exhale",...
"Hold",...
"Masker Noise ITD 50 Targ Right",...
"Masker Speech ITD 500 Targ Right",...
"Masker Speech ILD 10 Targ Left",...
"Masker Speech ITD 50 Targ Right",...
"Masker Noise ITD 500 Targ Right",...
"Masker Speech ITD 50 Targ Left",...
"Masker Noise ITD 500 Targ Left",...
"Masker Speech ITD 500 Targ Left",...
"Masker Speech ILD 10 Targ Right",...
"Masker Noise ILD 10 Targ Left",...
"Masker Noise ILD 10 Targ Right",...
"Masker Noise ITD 50 Targ Left"];

target_left_conditions = [6,9,10,11,13,15];
target_right_conditions = [4,5,7,8,12,14];
masker_speech_conditions = [5,6,7,9,11,12];
masker_noise_conditions = [4,8,10,13,14,15];
itd_50_conditions = [4,7,9,15];
itd_500_conditions = [5,8,10,11];
ild_10_conditions = [6,12,13,14];

%% Channel Key
% 1: S1_D1Hbo
% 2: S1_D1Hbr
% 3: S1_D2Hbo
% 4: S1_D2Hbr
% 5: S2_D1Hbo
% 6: S2_D1Hbr
% 7: S2_D2Hbo
% 8: S2_D2Hbr
% 9: S3_D3Hbo
% 10: S3_D3Hbr
% 11: S4_D3Hbo
% 12: S4_D3Hbr
% 13: S5_D4Hbo
% 14: S5_D4Hbr
% 15: S5_D5Hbo
% 16: S5_D5Hbr
% 17: S6_D4Hbo
% 18: S6_D4Hbr
% 19: S6_D5Hbo
% 20: S6_D5Hbr
% 21: S7_D6Hbo
% 22: S7_D6Hbr
% 23: S7_D7Hbo
% 24: S7_D7Hbr
% 25: S8_D6Hbo
% 26: S8_D6Hbr
% 27: S8_D7Hbo
% 28: S8_D7Hbr
channels = ["S1_D1Hbo",...
"S1_D1Hbr",...
"S1_D2Hbo",...
"S1_D2Hbr",...
"S2_D1Hbo",...
"S2_D1Hbr",...
"S2_D2Hbo",...
"S2_D2Hbr",...
"S3_D3Hbo",...
"S3_D3Hbr",...
"S4_D3Hbo",...
"S4_D3Hbr",...
"S5_D4Hbo",...
"S5_D4Hbr",...
"S5_D5Hbo",...
"S5_D5Hbr",...
"S6_D4Hbo",...
"S6_D4Hbr",...
"S6_D5Hbo",...
"S6_D5Hbr",...
"S7_D6Hbo",...
"S7_D6Hbr",...
"S7_D7Hbo",...
"S7_D7Hbr",...
"S8_D6Hbo",...
"S8_D6Hbr",...
"S8_D7Hbo",...
"S8_D7Hbr"];

hbo_channels = 1:2:27;
hbr_channels = 2:2:28;
left_stg_channels = 21:28; % S7_D6, S7_D7, S8_D6, S8_D7
right_stg_channels = 13:20; % S5_D4, S5_D5, S6_D4, S6_D5
left_dlpfc_channels = 1:8; % S1_D1, S1_D2, S2_D1, S2_D2
right_dlpfc_channels = 9:12; % S3_D3, S4_D3

%% Build Block Averages matrix
% Compile a numsubjects x numconditions x numchannels x num epochs X time matrix of
% baselined block average data
numconditions = length(conditions);
numchannels = length(channels);
block_average_data = [];
for isubject = 1:size(curr_subject_ID,2)
    curr_epochs = readtable(curr_subject_ID(isubject) + " block averages.csv");
    for icondition = 4:numconditions % excluding breath conditions
        epochs_this_condition = curr_epochs(string(curr_epochs.condition) == conditions(icondition),:);
        for ichannel = 1:numchannels
            epochs_this_channel = epochs_this_condition(:,[2,4,ichannel + 4]); % time, epoch number, and channel values
            unique_epochs = unique(epochs_this_channel.epoch);
            for iepoch = 1:length(unique_epochs) % for each epoch in this channel and condition
                % For each epoch, put the whole block average in the structure
                block_average_data(isubject,icondition,ichannel,iepoch,:) = table2array(epochs_this_channel(epochs_this_channel.epoch == unique_epochs(iepoch),3));
            end
        end
    end
end


%% Build GLM data matrix
% Compile a numsubjects x numconditions x numchannels X numepochs matrix of GLM beta
% values


%% Noise Masker vs. Speech Masker 
% Left STG, Right STG, Left DLPFC, Right DLPFC

% Block Average
lower_ylim = -0.2;
upper_ylim = 0.3;
figure;
subplot(2,2,1) % Left STG
hold on
time = linspace(-5,25,size(block_average_data,5));
mean_noise_masker_hbo = squeeze(mean(block_average_data(:,masker_noise_conditions,intersect(hbo_channels,left_stg_channels),:,:),[1,2,3,4]));
mean_noise_masker_hbr = squeeze(mean(block_average_data(:,masker_noise_conditions,intersect(hbr_channels,left_stg_channels),:,:),[1,2,3,4]));
mean_speech_masker_hbo = squeeze(mean(block_average_data(:,masker_speech_conditions,intersect(hbo_channels,left_stg_channels),:,:),[1,2,3,4]));
mean_speech_masker_hbr = squeeze(mean(block_average_data(:,masker_speech_conditions,intersect(hbr_channels,left_stg_channels),:,:),[1,2,3,4]));
SEM_noise_masker_hbo = squeeze(std(mean(block_average_data(:,masker_noise_conditions,intersect(hbo_channels,left_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_noise_masker_hbr = squeeze(std(mean(block_average_data(:,masker_noise_conditions,intersect(hbr_channels,left_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_speech_masker_hbo =squeeze(std(mean(block_average_data(:,masker_speech_conditions,intersect(hbo_channels,left_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_speech_masker_hbr = squeeze(std(mean(block_average_data(:,masker_speech_conditions,intersect(hbr_channels,left_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));

p_hbo_noise = shadedErrorBar(time,mean_noise_masker_hbo,SEM_noise_masker_hbo,'lineProps','-r');
p_hbr_noise = shadedErrorBar(time,mean_noise_masker_hbr,SEM_noise_masker_hbr,'lineProps','--r'); 
p_hbo_speech = shadedErrorBar(time,mean_speech_masker_hbo,SEM_speech_masker_hbo,'lineProps','-g');
p_hbr_speech = shadedErrorBar(time,mean_speech_masker_hbr,SEM_speech_masker_hbr,'lineProps','--g'); 

xlim([time(1),time(end)])
ylim([lower_ylim upper_ylim])
ylabel('\DeltaHb')
xlabel('Time from stimulus onset (s)')
title('Left STG')

subplot(2,2,2) % Right STG
hold on
time = linspace(-5,25,size(block_average_data,5));
mean_noise_masker_hbo = squeeze(mean(block_average_data(:,masker_noise_conditions,intersect(hbo_channels,right_stg_channels),:,:),[1,2,3,4]));
mean_noise_masker_hbr = squeeze(mean(block_average_data(:,masker_noise_conditions,intersect(hbr_channels,right_stg_channels),:,:),[1,2,3,4]));
mean_speech_masker_hbo = squeeze(mean(block_average_data(:,masker_speech_conditions,intersect(hbo_channels,right_stg_channels),:,:),[1,2,3,4]));
mean_speech_masker_hbr = squeeze(mean(block_average_data(:,masker_speech_conditions,intersect(hbr_channels,right_stg_channels),:,:),[1,2,3,4]));
SEM_noise_masker_hbo = squeeze(std(mean(block_average_data(:,masker_noise_conditions,intersect(hbo_channels,right_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_noise_masker_hbr = squeeze(std(mean(block_average_data(:,masker_noise_conditions,intersect(hbr_channels,right_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_speech_masker_hbo =squeeze(std(mean(block_average_data(:,masker_speech_conditions,intersect(hbo_channels,right_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_speech_masker_hbr = squeeze(std(mean(block_average_data(:,masker_speech_conditions,intersect(hbr_channels,right_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
p_hbo_noise = shadedErrorBar(time,mean_noise_masker_hbo,SEM_noise_masker_hbo,'lineProps','-r');
p_hbr_noise = shadedErrorBar(time,mean_noise_masker_hbr,SEM_noise_masker_hbr,'lineProps','--r'); 
p_hbo_speech = shadedErrorBar(time,mean_speech_masker_hbo,SEM_speech_masker_hbo,'lineProps','-g');
p_hbr_speech = shadedErrorBar(time,mean_speech_masker_hbr,SEM_speech_masker_hbr,'lineProps','--g'); 

xlim([time(1),time(end)])
ylim([lower_ylim upper_ylim])
ylabel('\DeltaHb')
xlabel('Time from stimulus onset (s)')
title('Right STG')

subplot(2,2,3) % Left DLPFC
hold on
time = linspace(-5,25,size(block_average_data,5));
mean_noise_masker_hbo = squeeze(mean(block_average_data(:,masker_noise_conditions,intersect(hbo_channels,left_dlpfc_channels),:,:),[1,2,3,4]));
mean_noise_masker_hbr = squeeze(mean(block_average_data(:,masker_noise_conditions,intersect(hbr_channels,left_dlpfc_channels),:,:),[1,2,3,4]));
mean_speech_masker_hbo = squeeze(mean(block_average_data(:,masker_speech_conditions,intersect(hbo_channels,left_dlpfc_channels),:,:),[1,2,3,4]));
mean_speech_masker_hbr = squeeze(mean(block_average_data(:,masker_speech_conditions,intersect(hbr_channels,left_dlpfc_channels),:,:),[1,2,3,4]));
SEM_noise_masker_hbo = squeeze(std(mean(block_average_data(:,masker_noise_conditions,intersect(hbo_channels,left_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_noise_masker_hbr = squeeze(std(mean(block_average_data(:,masker_noise_conditions,intersect(hbr_channels,left_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_speech_masker_hbo =squeeze(std(mean(block_average_data(:,masker_speech_conditions,intersect(hbo_channels,left_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_speech_masker_hbr = squeeze(std(mean(block_average_data(:,masker_speech_conditions,intersect(hbr_channels,left_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
p_hbo_noise = shadedErrorBar(time,mean_noise_masker_hbo,SEM_noise_masker_hbo,'lineProps','-r');
p_hbr_noise = shadedErrorBar(time,mean_noise_masker_hbr,SEM_noise_masker_hbr,'lineProps','--r'); 
p_hbo_speech = shadedErrorBar(time,mean_speech_masker_hbo,SEM_speech_masker_hbo,'lineProps','-g');
p_hbr_speech = shadedErrorBar(time,mean_speech_masker_hbr,SEM_speech_masker_hbr,'lineProps','--g'); 

xlim([time(1),time(end)])
ylim([lower_ylim upper_ylim])
ylabel('\DeltaHb')
xlabel('Time from stimulus onset (s)')
title('Left DLPFC')

subplot(2,2,4) % Right DLPFC
hold on
time = linspace(-5,25,size(block_average_data,5));
mean_noise_masker_hbo = squeeze(mean(block_average_data(:,masker_noise_conditions,intersect(hbo_channels,right_dlpfc_channels),:,:),[1,2,3,4]));
mean_noise_masker_hbr = squeeze(mean(block_average_data(:,masker_noise_conditions,intersect(hbr_channels,right_dlpfc_channels),:,:),[1,2,3,4]));
mean_speech_masker_hbo = squeeze(mean(block_average_data(:,masker_speech_conditions,intersect(hbo_channels,right_dlpfc_channels),:,:),[1,2,3,4]));
mean_speech_masker_hbr = squeeze(mean(block_average_data(:,masker_speech_conditions,intersect(hbr_channels,right_dlpfc_channels),:,:),[1,2,3,4]));
SEM_noise_masker_hbo = squeeze(std(mean(block_average_data(:,masker_noise_conditions,intersect(hbo_channels,right_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_noise_masker_hbr = squeeze(std(mean(block_average_data(:,masker_noise_conditions,intersect(hbr_channels,right_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_speech_masker_hbo =squeeze(std(mean(block_average_data(:,masker_speech_conditions,intersect(hbo_channels,right_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_speech_masker_hbr = squeeze(std(mean(block_average_data(:,masker_speech_conditions,intersect(hbr_channels,right_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
p_hbo_noise = shadedErrorBar(time,mean_noise_masker_hbo,SEM_noise_masker_hbo,'lineProps','-r');
p_hbr_noise = shadedErrorBar(time,mean_noise_masker_hbr,SEM_noise_masker_hbr,'lineProps','--r'); 
p_hbo_speech = shadedErrorBar(time,mean_speech_masker_hbo,SEM_speech_masker_hbo,'lineProps','-g');
p_hbr_speech = shadedErrorBar(time,mean_speech_masker_hbr,SEM_speech_masker_hbr,'lineProps','--g'); 

xlim([time(1),time(end)])
ylim([lower_ylim upper_ylim])
ylabel('\DeltaHb')
xlabel('Time from stimulus onset (s)')
title('Right DLPFC')

sgtitle('Speech vs. Noise') 
legend({'HbO Noise','HbR Noise','HbO Speech','HbR Speech'})


%% Target Left vs. Target Right
% Left STG, Right STG, Left DLPFC, Right DLPFC

% Block Average
figure;
subplot(2,2,1) % Left STG
hold on
time = linspace(-5,25,size(block_average_data,5));
mean_target_left_hbo = squeeze(mean(block_average_data(:,target_left_conditions,intersect(hbo_channels,left_stg_channels),:,:),[1,2,3,4]));
mean_target_left_hbr = squeeze(mean(block_average_data(:,target_left_conditions,intersect(hbr_channels,left_stg_channels),:,:),[1,2,3,4]));
mean_target_right_hbo = squeeze(mean(block_average_data(:,target_right_conditions,intersect(hbo_channels,left_stg_channels),:,:),[1,2,3,4]));
mean_target_right_hbr = squeeze(mean(block_average_data(:,target_right_conditions,intersect(hbr_channels,left_stg_channels),:,:),[1,2,3,4]));
SEM_target_left_hbo = squeeze(std(mean(block_average_data(:,target_left_conditions,intersect(hbo_channels,left_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_target_left_hbr = squeeze(std(mean(block_average_data(:,target_left_conditions,intersect(hbr_channels,left_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_target_right_hbo =squeeze(std(mean(block_average_data(:,target_right_conditions,intersect(hbo_channels,left_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_target_right_hbr = squeeze(std(mean(block_average_data(:,target_right_conditions,intersect(hbr_channels,left_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
p_hbo_left = shadedErrorBar(time,mean_target_left_hbo,SEM_target_left_hbo,'lineProps','-r');
p_hbr_left = shadedErrorBar(time,mean_target_left_hbr,SEM_target_left_hbr,'lineProps','--r'); 
p_hbo_right = shadedErrorBar(time,mean_target_right_hbo,SEM_target_right_hbo,'lineProps','-g');
p_hbr_right = shadedErrorBar(time,mean_target_right_hbr,SEM_target_right_hbr,'lineProps','--g'); 

xlim([time(1),time(end)])
ylim([lower_ylim upper_ylim])
ylabel('\DeltaHb')
xlabel('Time from stimulus onset (s)')
title('Left STG')

subplot(2,2,2) % Right STG
hold on
time = linspace(-5,25,size(block_average_data,5));
mean_target_left_hbo = squeeze(mean(block_average_data(:,target_left_conditions,intersect(hbo_channels,right_stg_channels),:,:),[1,2,3,4]));
mean_target_left_hbr = squeeze(mean(block_average_data(:,target_left_conditions,intersect(hbr_channels,right_stg_channels),:,:),[1,2,3,4]));
mean_target_right_hbo = squeeze(mean(block_average_data(:,target_right_conditions,intersect(hbo_channels,right_stg_channels),:,:),[1,2,3,4]));
mean_target_right_hbr = squeeze(mean(block_average_data(:,target_right_conditions,intersect(hbr_channels,right_stg_channels),:,:),[1,2,3,4]));
SEM_target_left_hbo = squeeze(std(mean(block_average_data(:,target_left_conditions,intersect(hbo_channels,right_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_target_left_hbr = squeeze(std(mean(block_average_data(:,target_left_conditions,intersect(hbr_channels,right_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_target_right_hbo =squeeze(std(mean(block_average_data(:,target_right_conditions,intersect(hbo_channels,right_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_target_right_hbr = squeeze(std(mean(block_average_data(:,target_right_conditions,intersect(hbr_channels,right_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
p_hbo_left = shadedErrorBar(time,mean_target_left_hbo,SEM_target_left_hbo,'lineProps','-r');
p_hbr_left = shadedErrorBar(time,mean_target_left_hbr,SEM_target_left_hbr,'lineProps','--r'); 
p_hbo_right = shadedErrorBar(time,mean_target_right_hbo,SEM_target_right_hbo,'lineProps','-g');
p_hbr_right = shadedErrorBar(time,mean_target_right_hbr,SEM_target_right_hbr,'lineProps','--g'); 

xlim([time(1),time(end)])
ylim([lower_ylim upper_ylim])
ylabel('\DeltaHb')
xlabel('Time from stimulus onset (s)')
title('Right STG')


subplot(2,2,3) % Left DLPFC
hold on
time = linspace(-5,25,size(block_average_data,5));
mean_target_left_hbo = squeeze(mean(block_average_data(:,target_left_conditions,intersect(hbo_channels,left_dlpfc_channels),:,:),[1,2,3,4]));
mean_target_left_hbr = squeeze(mean(block_average_data(:,target_left_conditions,intersect(hbr_channels,left_dlpfc_channels),:,:),[1,2,3,4]));
mean_target_right_hbo = squeeze(mean(block_average_data(:,target_right_conditions,intersect(hbo_channels,left_dlpfc_channels),:,:),[1,2,3,4]));
mean_target_right_hbr = squeeze(mean(block_average_data(:,target_right_conditions,intersect(hbr_channels,left_dlpfc_channels),:,:),[1,2,3,4]));
SEM_target_left_hbo = squeeze(std(mean(block_average_data(:,target_left_conditions,intersect(hbo_channels,left_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_target_left_hbr = squeeze(std(mean(block_average_data(:,target_left_conditions,intersect(hbr_channels,left_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_target_right_hbo =squeeze(std(mean(block_average_data(:,target_right_conditions,intersect(hbo_channels,left_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_target_right_hbr = squeeze(std(mean(block_average_data(:,target_right_conditions,intersect(hbr_channels,left_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
p_hbo_left = shadedErrorBar(time,mean_target_left_hbo,SEM_target_left_hbo,'lineProps','-r');
p_hbr_left = shadedErrorBar(time,mean_target_left_hbr,SEM_target_left_hbr,'lineProps','--r'); 
p_hbo_right = shadedErrorBar(time,mean_target_right_hbo,SEM_target_right_hbo,'lineProps','-g');
p_hbr_right = shadedErrorBar(time,mean_target_right_hbr,SEM_target_right_hbr,'lineProps','--g'); 

xlim([time(1),time(end)])
ylim([lower_ylim upper_ylim])
ylabel('\DeltaHb')
xlabel('Time from stimulus onset (s)')
title('Left DLPFC')

subplot(2,2,4) % Right DLPFC
hold on
time = linspace(-5,25,size(block_average_data,5));
mean_target_left_hbo = squeeze(mean(block_average_data(:,target_left_conditions,intersect(hbo_channels,right_dlpfc_channels),:,:),[1,2,3,4]));
mean_target_left_hbr = squeeze(mean(block_average_data(:,target_left_conditions,intersect(hbr_channels,right_dlpfc_channels),:,:),[1,2,3,4]));
mean_target_right_hbo = squeeze(mean(block_average_data(:,target_right_conditions,intersect(hbo_channels,right_dlpfc_channels),:,:),[1,2,3,4]));
mean_target_right_hbr = squeeze(mean(block_average_data(:,target_right_conditions,intersect(hbr_channels,right_dlpfc_channels),:,:),[1,2,3,4]));
SEM_target_left_hbo = squeeze(std(mean(block_average_data(:,target_left_conditions,intersect(hbo_channels,right_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_target_left_hbr = squeeze(std(mean(block_average_data(:,target_left_conditions,intersect(hbr_channels,right_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_target_right_hbo =squeeze(std(mean(block_average_data(:,target_right_conditions,intersect(hbo_channels,right_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_target_right_hbr = squeeze(std(mean(block_average_data(:,target_right_conditions,intersect(hbr_channels,right_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
p_hbo_left = shadedErrorBar(time,mean_target_left_hbo,SEM_target_left_hbo,'lineProps','-r');
p_hbr_left = shadedErrorBar(time,mean_target_left_hbr,SEM_target_left_hbr,'lineProps','--r'); 
p_hbo_right = shadedErrorBar(time,mean_target_right_hbo,SEM_target_right_hbo,'lineProps','-g');
p_hbr_right = shadedErrorBar(time,mean_target_right_hbr,SEM_target_right_hbr,'lineProps','--g'); 


xlim([time(1),time(end)])
ylim([lower_ylim upper_ylim])
ylabel('\DeltaHb')
xlabel('Time from stimulus onset (s)')
title('Right DLPFC')

sgtitle('Target Left vs. Target Right') 
legend({'HbO Target Left','HbR Target Left','HbO Target Right','HbR Target Right'})

%% ITD 50 vs. ITD 500 vs. ILD 10 
% Left STG, Right STG, Left DLPFC, Right DLPFC

% Block Average
figure;
subplot(2,2,1) % Left STG
hold on

mean_itd_small_hbo = squeeze(mean(block_average_data(:,itd_50_conditions,intersect(hbo_channels,left_stg_channels),:,:),[1,2,3,4]));
mean_itd_small_hbr = squeeze(mean(block_average_data(:,itd_50_conditions,intersect(hbr_channels,left_stg_channels),:,:),[1,2,3,4]));
mean_itd_large_hbo = squeeze(mean(block_average_data(:,itd_500_conditions,intersect(hbo_channels,left_stg_channels),:,:),[1,2,3,4]));
mean_itd_large_hbr = squeeze(mean(block_average_data(:,itd_500_conditions,intersect(hbr_channels,left_stg_channels),:,:),[1,2,3,4]));
mean_ild_hbo = squeeze(mean(block_average_data(:,ild_10_conditions,intersect(hbo_channels,left_stg_channels),:,:),[1,2,3,4]));
mean_ild_hbr = squeeze(mean(block_average_data(:,ild_10_conditions,intersect(hbr_channels,left_stg_channels),:,:),[1,2,3,4]));
SEM_itd_small_hbo = squeeze(std(mean(block_average_data(:,itd_50_conditions,intersect(hbo_channels,left_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_itd_small_hbr = squeeze(std(mean(block_average_data(:,itd_50_conditions,intersect(hbr_channels,left_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_itd_large_hbo =squeeze(std(mean(block_average_data(:,itd_500_conditions,intersect(hbo_channels,left_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_itd_large_hbr = squeeze(std(mean(block_average_data(:,itd_500_conditions,intersect(hbr_channels,left_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_ild_hbo =squeeze(std(mean(block_average_data(:,ild_10_conditions,intersect(hbo_channels,left_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_ild_hbr = squeeze(std(mean(block_average_data(:,ild_10_conditions,intersect(hbr_channels,left_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
p_hbo_itd_small = shadedErrorBar(time,mean_itd_small_hbo,SEM_itd_small_hbo,'lineProps','-r');
p_hbr_itd_small = shadedErrorBar(time,mean_itd_small_hbr,SEM_itd_small_hbr,'lineProps','--r'); 
p_hbo_itd_large = shadedErrorBar(time,mean_itd_large_hbo,SEM_itd_large_hbo,'lineProps','-g');
p_hbr_itd_large = shadedErrorBar(time,mean_itd_large_hbr,SEM_itd_large_hbr,'lineProps','--g'); 
p_hbo_ild = shadedErrorBar(time,mean_ild_hbo,SEM_ild_hbo,'lineProps','-b');
p_hbr_ild = shadedErrorBar(time,mean_ild_hbr,SEM_ild_hbr,'lineProps','--b'); 

xlim([time(1),time(end)])
ylim([lower_ylim upper_ylim])
ylabel('\DeltaHb')
xlabel('Time from stimulus onset (s)')
title('Left STG')

subplot(2,2,2) % Left STG
hold on

mean_itd_small_hbo = squeeze(mean(block_average_data(:,itd_50_conditions,intersect(hbo_channels,right_stg_channels),:,:),[1,2,3,4]));
mean_itd_small_hbr = squeeze(mean(block_average_data(:,itd_50_conditions,intersect(hbr_channels,right_stg_channels),:,:),[1,2,3,4]));
mean_itd_large_hbo = squeeze(mean(block_average_data(:,itd_500_conditions,intersect(hbo_channels,right_stg_channels),:,:),[1,2,3,4]));
mean_itd_large_hbr = squeeze(mean(block_average_data(:,itd_500_conditions,intersect(hbr_channels,right_stg_channels),:,:),[1,2,3,4]));
mean_ild_hbo = squeeze(mean(block_average_data(:,ild_10_conditions,intersect(hbo_channels,right_stg_channels),:,:),[1,2,3,4]));
mean_ild_hbr = squeeze(mean(block_average_data(:,ild_10_conditions,intersect(hbr_channels,right_stg_channels),:,:),[1,2,3,4]));
SEM_itd_small_hbo = squeeze(std(mean(block_average_data(:,itd_50_conditions,intersect(hbo_channels,right_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_itd_small_hbr = squeeze(std(mean(block_average_data(:,itd_50_conditions,intersect(hbr_channels,right_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_itd_large_hbo =squeeze(std(mean(block_average_data(:,itd_500_conditions,intersect(hbo_channels,right_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_itd_large_hbr = squeeze(std(mean(block_average_data(:,itd_500_conditions,intersect(hbr_channels,right_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_ild_hbo =squeeze(std(mean(block_average_data(:,ild_10_conditions,intersect(hbo_channels,right_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_ild_hbr = squeeze(std(mean(block_average_data(:,ild_10_conditions,intersect(hbr_channels,right_stg_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
p_hbo_itd_small = shadedErrorBar(time,mean_itd_small_hbo,SEM_itd_small_hbo,'lineProps','-r');
p_hbr_itd_small = shadedErrorBar(time,mean_itd_small_hbr,SEM_itd_small_hbr,'lineProps','--r'); 
p_hbo_itd_large = shadedErrorBar(time,mean_itd_large_hbo,SEM_itd_large_hbo,'lineProps','-g');
p_hbr_itd_large = shadedErrorBar(time,mean_itd_large_hbr,SEM_itd_large_hbr,'lineProps','--g'); 
p_hbo_ild = shadedErrorBar(time,mean_ild_hbo,SEM_ild_hbo,'lineProps','-b');
p_hbr_ild = shadedErrorBar(time,mean_ild_hbr,SEM_ild_hbr,'lineProps','--b'); 

xlim([time(1),time(end)])
ylim([lower_ylim upper_ylim])
ylabel('\DeltaHb')
xlabel('Time from stimulus onset (s)')
title('Right STG')

subplot(2,2,3) % Left STG
hold on

mean_itd_small_hbo = squeeze(mean(block_average_data(:,itd_50_conditions,intersect(hbo_channels,left_dlpfc_channels),:,:),[1,2,3,4]));
mean_itd_small_hbr = squeeze(mean(block_average_data(:,itd_50_conditions,intersect(hbr_channels,left_dlpfc_channels),:,:),[1,2,3,4]));
mean_itd_large_hbo = squeeze(mean(block_average_data(:,itd_500_conditions,intersect(hbo_channels,left_dlpfc_channels),:,:),[1,2,3,4]));
mean_itd_large_hbr = squeeze(mean(block_average_data(:,itd_500_conditions,intersect(hbr_channels,left_dlpfc_channels),:,:),[1,2,3,4]));
mean_ild_hbo = squeeze(mean(block_average_data(:,ild_10_conditions,intersect(hbo_channels,left_dlpfc_channels),:,:),[1,2,3,4]));
mean_ild_hbr = squeeze(mean(block_average_data(:,ild_10_conditions,intersect(hbr_channels,left_dlpfc_channels),:,:),[1,2,3,4]));
SEM_itd_small_hbo = squeeze(std(mean(block_average_data(:,itd_50_conditions,intersect(hbo_channels,left_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_itd_small_hbr = squeeze(std(mean(block_average_data(:,itd_50_conditions,intersect(hbr_channels,left_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_itd_large_hbo =squeeze(std(mean(block_average_data(:,itd_500_conditions,intersect(hbo_channels,left_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_itd_large_hbr = squeeze(std(mean(block_average_data(:,itd_500_conditions,intersect(hbr_channels,left_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_ild_hbo =squeeze(std(mean(block_average_data(:,ild_10_conditions,intersect(hbo_channels,left_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_ild_hbr = squeeze(std(mean(block_average_data(:,ild_10_conditions,intersect(hbr_channels,left_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
p_hbo_itd_small = shadedErrorBar(time,mean_itd_small_hbo,SEM_itd_small_hbo,'lineProps','-r');
p_hbr_itd_small = shadedErrorBar(time,mean_itd_small_hbr,SEM_itd_small_hbr,'lineProps','--r'); 
p_hbo_itd_large = shadedErrorBar(time,mean_itd_large_hbo,SEM_itd_large_hbo,'lineProps','-g');
p_hbr_itd_large = shadedErrorBar(time,mean_itd_large_hbr,SEM_itd_large_hbr,'lineProps','--g'); 
p_hbo_ild = shadedErrorBar(time,mean_ild_hbo,SEM_ild_hbo,'lineProps','-b');
p_hbr_ild = shadedErrorBar(time,mean_ild_hbr,SEM_ild_hbr,'lineProps','--b'); 

xlim([time(1),time(end)])
ylim([lower_ylim upper_ylim])
ylabel('\DeltaHb')
xlabel('Time from stimulus onset (s)')
title('Left DLPFC')

subplot(2,2,4) % Left STG
hold on

mean_itd_small_hbo = squeeze(mean(block_average_data(:,itd_50_conditions,intersect(hbo_channels,right_dlpfc_channels),:,:),[1,2,3,4]));
mean_itd_small_hbr = squeeze(mean(block_average_data(:,itd_50_conditions,intersect(hbr_channels,right_dlpfc_channels),:,:),[1,2,3,4]));
mean_itd_large_hbo = squeeze(mean(block_average_data(:,itd_500_conditions,intersect(hbo_channels,right_dlpfc_channels),:,:),[1,2,3,4]));
mean_itd_large_hbr = squeeze(mean(block_average_data(:,itd_500_conditions,intersect(hbr_channels,right_dlpfc_channels),:,:),[1,2,3,4]));
mean_ild_hbo = squeeze(mean(block_average_data(:,ild_10_conditions,intersect(hbo_channels,right_dlpfc_channels),:,:),[1,2,3,4]));
mean_ild_hbr = squeeze(mean(block_average_data(:,ild_10_conditions,intersect(hbr_channels,right_dlpfc_channels),:,:),[1,2,3,4]));
SEM_itd_small_hbo = squeeze(std(mean(block_average_data(:,itd_50_conditions,intersect(hbo_channels,right_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_itd_small_hbr = squeeze(std(mean(block_average_data(:,itd_50_conditions,intersect(hbr_channels,right_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_itd_large_hbo =squeeze(std(mean(block_average_data(:,itd_500_conditions,intersect(hbo_channels,right_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_itd_large_hbr = squeeze(std(mean(block_average_data(:,itd_500_conditions,intersect(hbr_channels,right_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_ild_hbo =squeeze(std(mean(block_average_data(:,ild_10_conditions,intersect(hbo_channels,right_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
SEM_ild_hbr = squeeze(std(mean(block_average_data(:,ild_10_conditions,intersect(hbr_channels,right_dlpfc_channels),:,:),[2,4]),[],[1,2,3,4]))./(sqrt(length(curr_subject_ID)));
p_hbo_itd_small = shadedErrorBar(time,mean_itd_small_hbo,SEM_itd_small_hbo,'lineProps','-r');
p_hbr_itd_small = shadedErrorBar(time,mean_itd_small_hbr,SEM_itd_small_hbr,'lineProps','--r'); 
p_hbo_itd_large = shadedErrorBar(time,mean_itd_large_hbo,SEM_itd_large_hbo,'lineProps','-g');
p_hbr_itd_large = shadedErrorBar(time,mean_itd_large_hbr,SEM_itd_large_hbr,'lineProps','--g'); 
p_hbo_ild = shadedErrorBar(time,mean_ild_hbo,SEM_ild_hbo,'lineProps','-b');
p_hbr_ild = shadedErrorBar(time,mean_ild_hbr,SEM_ild_hbr,'lineProps','--b'); 

xlim([time(1),time(end)])
ylim([lower_ylim upper_ylim])
ylabel('\DeltaHb')
xlabel('Time from stimulus onset (s)')
title('Right DLPFC')


sgtitle('Comparison Across Spatial Cues') 
legend({'HbO ITD Small','HbR ITD Small','HbO ITD Large','HbR ITD Large','HbO ILD','HbR ILD'})



%% Breakdown by subject
upper_ylim = 0.4;
lower_ylim = -0.2;
possible_colors = ["r","g","b","y"];
for isubject = 1:size(block_average_data,1)
    figure;
    all_left_stg_blocks_hbo = squeeze(block_average_data(isubject,4:end,intersect(hbo_channels,left_stg_channels),:,:));
    subplot(2,2,1)
    hold on
    for ichannel = 1:size(all_left_stg_blocks_hbo,2)
            curr_color = possible_colors(ichannel);
                shadedErrorBar(time,squeeze(mean(all_left_stg_blocks_hbo(:,ichannel,:,:),[1,3])),squeeze(std(all_left_stg_blocks_hbo(:,ichannel,:,:),[],[1,3]))./sqrt(size(all_left_stg_blocks_hbo,3)),'lineProps',char(curr_color))
        end
    xlim([time(1),time(end)])
    ylim([lower_ylim upper_ylim])
    ylabel('\DeltaHbO')
    xlabel('Time from stimulus onset (s)')
    title('Left STG')

    all_right_stg_blocks_hbo = squeeze(block_average_data(isubject,4:end,intersect(hbo_channels,right_stg_channels),:,:));
    subplot(2,2,2)
    hold on
    for ichannel = 1:size(all_right_stg_blocks_hbo,2)
            curr_color = possible_colors(ichannel);
                shadedErrorBar(time,squeeze(mean(all_right_stg_blocks_hbo(:,ichannel,:,:),[1,3])),squeeze(std(all_right_stg_blocks_hbo(:,ichannel,:,:),[],[1,3]))./sqrt(size(all_right_stg_blocks_hbo,3)),'lineProps',char(curr_color))
        end
    xlim([time(1),time(end)])
    ylim([lower_ylim upper_ylim])
    ylabel('\DeltaHbO')
    xlabel('Time from stimulus onset (s)')
    title('Right STG')

    all_left_dlpfc_blocks_hbo = squeeze(block_average_data(isubject,4:end,intersect(hbo_channels,left_dlpfc_channels),:,:));
    subplot(2,2,3)
    hold on
    for ichannel = 1:size(all_left_dlpfc_blocks_hbo,2)
            curr_color = possible_colors(ichannel);
                shadedErrorBar(time,squeeze(mean(all_left_dlpfc_blocks_hbo(:,ichannel,:,:),[1,3])),squeeze(std(all_left_dlpfc_blocks_hbo(:,ichannel,:,:),[],[1,3]))./sqrt(size(all_left_dlpfc_blocks_hbo,3)),'lineProps',char(curr_color))
    end
    xlim([time(1),time(end)])
    ylim([lower_ylim upper_ylim])
    ylabel('\DeltaHbO')
    xlabel('Time from stimulus onset (s)')
    title('Left DLPFC')

    all_right_dlpfc_blocks_hbo = squeeze(block_average_data(isubject,4:end,intersect(hbo_channels,right_dlpfc_channels),:,:));
    subplot(2,2,4)
    hold on
        for ichannel = 1:size(all_right_dlpfc_blocks_hbo,2)
            curr_color = possible_colors(ichannel);
                shadedErrorBar(time,squeeze(mean(all_right_dlpfc_blocks_hbo(:,ichannel,:,:),[1,3])),squeeze(std(all_right_dlpfc_blocks_hbo(:,ichannel,:,:),[],[1,3]))./sqrt(size(all_right_dlpfc_blocks_hbo,3)),'lineProps',char(curr_color))
        end
    xlim([time(1),time(end)])
    ylim([lower_ylim upper_ylim])
    ylabel('\DeltaHbO')
    xlabel('Time from stimulus onset (s)')
    title('Right DLPFC')

    sgtitle(curr_subject_ID(isubject))
end
