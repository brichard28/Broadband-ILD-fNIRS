%% BNR_master_script_for_computing_node

% This script calls other analysis functions in order to run a script on
% the cluster node for the Oscillations project.
%% Define directory
cd /lab_data/barblab/Ben/Oscillations
addpath('/lab_data/barblab/Ben/Oscillations')
addpath('/lab_data/barblab/Ben/Oscillations/Fieldtrip-lite170623')
addpath('/lab_data/barblab/Ben/Oscillations/BNR Alpha Analyses')
addpath('/lab_data/barblab/Ben/Oscillations/BNR Final DATA Files')
addpath('/lab_data/barblab/Ben/Oscillations/eeglab2021.1')

%BNR_group_analysis_processing_bytrial;

curr_time = clock;
results_filename = join(['Oscillations_Results_BNR_',join(split(string(curr_time(1:end-1)),','),'_')]);
results_filename = join([results_filename,'.txt'],'');
fileID = fopen(results_filename,'w');

%% Using Only Day 1 Data
subject_tags = {'018XC','019XX','021XC','022XX','023XX','025MX','026MC','027MC',...
    '028MC','029MX','030MX','031MC',...
    '031XX','032MX','034MX','035MX','039XC','040MC',...
    '051XC','053XX','055XX','061XX','063XX','064XC','067XX',...
    '068XX','070XX','073XX','075XC','076XX','080XX','084MX',...
    '091MC','093XX','096MC','098MX',...
    '100MX','101XC','102XX',...
    '103MX','105MX','110XX',...
    '112MX','200XC'};

fwrite(fileID,['The subjects run in this analysis are:' newline]);
fwrite(fileID,join(string(subject_tags),','));

%% Find Best Alpha Frequency
best_alpha_freqs = [13,12.75,10.25,10.25,12.5,9.25,10.75,11.75,10.75,9.5,10.75,14,14,11,10.75,14,14,13.25,10.25,10.5,10.25,10.75,8,10.75,8,10,9.75,11.25,11.25,10.5,9.75,10.25,10.5,13.25,14,11,10.25,11.5,11.5,11,9,11,10,12];
% This will give the IPAF for each subject in variable best_alpha_freqs
fwrite(fileID,newline);
fwrite(fileID,[newline, 'The IPAF For Each Subject is:', newline]);
fwrite(fileID,join(string(best_alpha_freqs),','));

%% Find is_adhd vector from excel file of Subject Info
subject_info = readtable('SubjectInfo.xlsx','Sheet','Subject Info');
all_subject_IDs = subject_info.(1); % subject tags from the excel sheet
is_adhd = subject_info.(2); % 1s or 0s for ADHD
is_adhd = is_adhd(ismember(all_subject_IDs,subject_tags));
analysis_info = readtable('SubjectInfo.xlsx','Sheet','Analysis Info');

% %% Calculate Spectrogram by trial
% BNR_group_analysis_processing_bytrial;

%% Find alpha for every file, store it
LFC_channels = 1:11;
LPO_channels = 16:27;
RFC_channels = 32 + [2:4,7:14];
RPO_channels = 32 + [21:32];
APO_channels = [LPO_channels,28:32,RPO_channels];
AFC_channels = [LFC_channels,33,32+5,32+6,32+15,RFC_channels];


%subject_tags = subject_tags(is_adhd == 1);
channel_select = APO_channels;

parfor isubject = 1:length(subject_tags)
    spectrogram = load(string(append(subject_tags(isubject),'_spectrogram_bytrial a la Shim.mat')));
    power = spectrogram.power;
    freqs = spectrogram.freqs;
    times = spectrogram.times;
    curr_subj_tag = split(string(subject_tags(isubject)),'');
    curr_subj_tag = curr_subj_tag(2:6)';
    scoreanderp = load(join(string(['AUD_ASA',curr_subj_tag(1:3),'_',curr_subj_tag(4:5),'_1-50_DATA.mat']),''));
    SCORE = scoreanderp.SCORE;
    ERP = scoreanderp.ERP;
    working_power = power;
    %all_powers(isubject,:,:,:,:) = squeeze(working_power(:,:,:,:));
    current_alpha_power = squeeze(mean(working_power(:,:,find(freqs==best_alpha_freqs(isubject)-2):find(freqs==best_alpha_freqs(isubject)+2),:),[3])); % take average over frequency
    curr_size = size(current_alpha_power);
    disp(curr_size(1));
    if curr_size(1) < 240
        current_alpha_power(curr_size(1):240,:,:) = nan;
    end
    all_alpha_powers(isubject,:,:,:) = current_alpha_power; % this should be trial by time
    [~,timeindex1] = min(abs(times(:)+1600));
    [~,timeindex2] = min(abs(times(:)-0));
    all_alpha_powers_baselined(isubject,:,:,:) = current_alpha_power; %squeeze(((current_alpha_power - mean(current_alpha_power(:,:,1:timeindex1),'all')))./(std(current_alpha_power(:,:,1:timeindex1),[],'all')));
end
freqsandtimes = load(string(append(subject_tags(1),'_spectrogram.mat')),'freqs','times');
freqs = freqsandtimes.freqs;
times = freqsandtimes.times;

% save that fucker
alpha_filename = join(['/lab_data/barblab/Ben/Oscillations/Alpha_Powers_BNR_',join(split(string(curr_time(1:end-1)),','),'_')]);
save(alpha_filename,'all_alpha_powers_baselined','freqs','times','-v7.3')
% fwrite(fileID,newline);
% fwrite(fileID,join(['All alpha powers is stored in the file',string(alpha_filename)]));

%clear all_alpha_powers_baselined
%% Cluster based analysis of spectrograms with variables electrode, frequency, and time
% [stat,avg] = BNR_cluster_analysis_main();
% %save('/lab_data/barblab/Ben/Oscillations/BNR_cluster_analysis_results','stat','avg');


%% Correspond with behavior comparisons
%Variance in alpha signal within trial types (version 2)
% [SEM_allsubjects_overtime,SEM_allsubjects,SEM_F_N_overtime,SEM_F_N,SEM_F_S1_overtime,SEM_F_S1,SEM_F_S2_overtime,SEM_F_S2,...
%     SEM_B_N_overtime,SEM_B_N,SEM_B_S1_overtime,SEM_B_S1,SEM_B_S2_overtime,SEM_B_S2,average_percent_correct,SEM_FOCAL,SEM_BROAD] = BNR_alpha_correspond_with_behavior_v2(best_alpha_freqs,all_alpha_powers_baselined);
% save('/lab_data/barblab/Ben/Oscillations/OscillationsSEM_results.mat','SEM_allsubjects_overtime','SEM_allsubjects','SEM_F_N_overtime','SEM_F_N','SEM_F_S1_overtime','SEM_F_S1','SEM_F_S2_overtime','SEM_F_S2',...
%     'SEM_B_N_overtime','SEM_B_N','SEM_B_S1_overtime','SEM_B_S1','SEM_B_S2_overtime','SEM_B_S2','average_percent_correct','SEM_FOCAL','SEM_BROAD');

% 
% fwrite(fileID,newline)
% fwrite(fileID,'The v2 correspond with behavior plots are stored in pdfs end in name _compare_with_behavior.pdf')
% fwrite(fileID,[newline,'The r values (in order All trials, F_N, F_S1, F_S2,B_N, B_S1, B_S2, are:',newline])
% %fwrite(fileID,join(string(r),','))
% fwrite(fileID,[newline,'The p values (in order All trials, F_N, F_S1, F_S2,B_N, B_S1, B_S2, are:',newline])
% %fwrite(fileID,join(string(p),','))

%% Output to a pdf

fclose(fileID);
