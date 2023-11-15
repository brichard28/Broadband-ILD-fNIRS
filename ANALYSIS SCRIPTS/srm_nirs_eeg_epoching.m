%% srm_nirs_eeg_epoching

% Benjamin Richardson
% Script to epoch EEG data for SRM-NIRS-EEG-1

% 63999 = trigger for event
addpath('C:\Users\benri\Documents\eeglab2023.0');
addpath('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\prepro_epoched_data')
all_subID = ["NDARZD647HJ1","NDARVX375BR6",'NDARHN971WJ5','NDARHM932KNX','NDARGF569BF3','NDARBL382XK5'];
subID = 'NDARBL382XK5';
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename',[subID, '_ICAdone.set'],'filepath','C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\prepro_epoched_data\');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset(EEG);
% Adjust event latencies to match with audio downsample
fs = EEG.srate;
delay = fs/44100;
shifting_latencies = mat2cell(cell2mat({EEG.event.latency}') + (delay * fs), length(EEG.event), 1);
shifting_latencies = shifting_latencies{:};
for i = 1:numel(shifting_latencies)
    EEG.event(i).latency = shifting_latencies(i);
end
EEG = eeg_checkset(EEG);
% All epochs
EEG = pop_epoch( EEG, {'63999'}, [-3 12], 'newname', [subID, 'All epochs'], 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off');
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [],[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off');
EEG = eeg_checkset( EEG );
save(['C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\prepro_epoched_data\' subID 'all_epochs.mat'],"EEG")
