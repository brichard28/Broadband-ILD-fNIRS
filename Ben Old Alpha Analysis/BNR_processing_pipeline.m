%% BNR_processing_pipeline
% Author: Ben Richardson
% Date Created: 10/7/2021

% Scripts this code calls/steps it takes:
% BNR_preprocessing_step1: load in both EEG BDF files for the subject.
% Remove unwanted channels, downsample, and sort into visual and auditory.
% Output: un-epoched/no ICA just visual and just auditory EEG


% BNR_preprocessing_step2: Run for visual and auditory separately. BP
% filter the data, add channel locations, run ICA, epoch, and reject
% components and weird trials. Output: cleaned up and pruned EEG data,
% which components were rejected, which trials are marked as weird

% BNR_combine_with_behavior: Output: _DATA file which contains structures
% ERP {subjID, erp EEG data, time vector, time limits, fs,
% eyeBlinkRej (1 or 0), baselinecorrected (1 or 0), numBadTrialsRemoved,
% badTrials (logical vector), filter frequencies, and date analyzed}
% and SCORE {SubjID, Stimulant Status, Caffeine Status, ADHD Status, ADHD
% subtype, perblockPC (per block percent correct), attend (1=attend target,
% 2 = attend supertarget), answ (correct answer in order), resp (subject
% response in order), condition ( 4 = supertarget, 3 = no supertarget),
% delay (1 = 1 second, 1.5 = 1.5 seconds), hits (24 x 10 logical),
% hitspertrial (24 x 3 x 10 logical), trial type specific logical for
% marking trial type (e.g. F_N), start time (block start times), stop time
% (block end times), respTime (response time within trial), parameters
% structure (totalDuration of trial, onDuration, timestamp, numBlocks,
% numTrials), dateprocessed, numValidTrials}
% ** Reminder: There are 10 blocks of 24 trials**

% BNR_plot_erps: Input: channels to average over, output: erp plots

% BNR_plot_time_freq: Input: newtimef function parameters, output:
% spectrogram plots (used for preliminary alpha analysis).

files_to_load = string(input('Please enter the files you would like to analayze (form: file1.bdf,file2.bdf): ','s'));
files_to_load = strsplit(files_to_load, ',');
file_for_info = char(files_to_load(1));
subject_ID = file_for_info(5:7);
stimstatus = file_for_info(9);
caffstatus = file_for_info(10);
BPF = '1-50'; % Make sure matches filter_edge_freqs 
%% Set parameters and run step 1
Fs = 256;
filter_edge_freqs = [1 50];
BNR_preprocessing_step1

%% Setparameters and run step 2
epoch_time_limits_vis =  [-2.5 4.5];
epoch_time_limits_aud = [-2.5 4.5];
% preprocess visual
EEG = EEG_visual;
epoch_time_limits = epoch_time_limits_vis;
BNR_preprocessing_step2;
EEG_visual = EEG;
[ALLEEG, EEG_visual, CURRENTSET] = pop_newset(ALLEEG, EEG_visual,CURRENTSET,'setname',append(dataset_base_name,'Visual Preprocessed'),'gui','off');

disp('Press Enter to continue')
input('')
% preprocess auditory
EEG = EEG_auditory;
epoch_time_limits = epoch_time_limits_aud;
BNR_preprocessing_step2;
EEG_auditory = EEG; 
[ALLEEG, EEG_auditory, CURRENTSET] = pop_newset(ALLEEG, EEG_auditory, CURRENTSET,'setname',append(dataset_base_name,'Auditory Preprocessed'),'gui','off');

%% Save preprocessed EEG
[FILENAME, PATHNAME, FILTERINDEX] = uiputfile(strcat('PREPROCESSED_VIS_ASA',subject_ID,'_',stimstatus,caffstatus,'_',BPF,'_EEG.mat'));
save([PATHNAME,FILENAME],'EEG_visual');

[FILENAME, PATHNAME, FILTERINDEX] = uiputfile(strcat('PREPROCESSED_AUD_ASA',subject_ID,'_',stimstatus,caffstatus,'_',BPF,'_EEG.mat'));
save([PATHNAME,FILENAME],'EEG_auditory');

%% Combine each EEG structure with behavior to get a _DATA.mat file for auditory and visual separately
is_vis = 1;
is_aud = 0;
epoch_time_limits = epoch_time_limits_vis;
EEG = EEG_visual;
BNR_combine_with_behavior
is_vis = 0;
is_aud = 1;
epoch_time_limits = epoch_time_limits_aud;
EEG = EEG_auditory;
BNR_combine_with_behavior
