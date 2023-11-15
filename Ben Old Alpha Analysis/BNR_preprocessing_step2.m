%% BNR_preprocessing_step2
% BNR_preprocessing_step2: Run for visual and auditory separately. BP
% filter the data, add channel locations, run ICA, epoch, and reject
% components and weird trials. Output: cleaned up and pruned EEG data,
% which components were rejected, which trials are marked as weird


% Bandpass filter the data (1-20 Hz Windowed sinc Kaiser for now)
EEG = pop_firws(EEG, 'fcutoff', filter_edge_freqs, 'ftype', 'bandpass', 'wtype', 'kaiser', 'warg', 5.65326, 'forder', 1856, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);
% Save BP filtered data as a new dataset
dataset_BPF_name = append(dataset_base_name, ' BPF');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',dataset_BPF_name,'gui','off');
%EEG = eeg_checkset(EEG);

% THIS IS WHERE YOU WOULD CLEAN UP, but we skip for now

% Add channel location data (from chanlocs_64_3_eye_chan.locs)
EEG=pop_chanedit(EEG, 'load',{'C:\Users\benri\Documents\PhD Year 1\Jasmine EEG Data\chanlocs_64_3_eye_chan.locs','filetype','loc'});
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%EEG = eeg_checkset(EEG);

% Run ICA Analysis, save as new set
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1);
%[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
dataset_ICA_name = append(dataset_BPF_name, ' ICA');
%EEG = pop_editset(EEG, 'setname',dataset_ICA_name, 'icaweights', 'EEG10.icaweights', 'icasphere', 'EEG10.icasphere', 'icachansind', 'EEG10.icachansind'); % take the ICA weights and apply them to this set
%EEG = eeg_checkset(EEG);

% Relabel events?? So that some of the triggers which mean the same thing
% are fine 

% Epoch the data (~2.5 seconds before onset of stim, and 4.5 seconds after)
dataset_epoched_name = append(dataset_ICA_name,' epochs');
EEG = pop_epoch(EEG, {'3841'}, epoch_time_limits, 'newname', dataset_epoched_name, 'epochinfo', 'yes'); %65281 3841
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off');
%EEG = eeg_checkset(EEG);
% Pick and remove bad components:
% TIPS: use topography, spectrum, and component activations. Brain data is
% usually diffuse (volume conduction), and the brain and body work at lower
% frequencies

pop_selectcomps(EEG,[1:64]);
pop_eegplot( EEG, 1, 1, 1); % plot channel activations (scroll)
pop_eegplot( EEG, 0, 1, 1); % plot component activations (scroll)
pause
%pop_prop( EEG, 0, 1, NaN, {'freqrange',[1 50] });


% prompt for which components should be removed
components_to_remove = input('Input which components you would like to remove, separated by commas (ex: [1,2,7]): ','s');
close all
components_to_remove = str2num(components_to_remove);
% Transfer ICA pruned data to another dataset (remove rejected components)
EEG = pop_subcomp( EEG, components_to_remove, 0);
dataset_pruned_name = append(dataset_epoched_name, ' pruned');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname', dataset_pruned_name,'gui','off'); 


%% Go through and manually mark weird trials. Save rej_manual 
pop_eegplot( EEG, 1, 1, 0);
disp('Sort through and manually mark weird trials. These will be stored in EEG.reject.rejmanual. Then, press Enter to continue.');

input('')

