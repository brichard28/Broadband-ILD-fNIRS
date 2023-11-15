%% BNR_preprocessing_step1
% This function is Ben's first attempt at building an EEG processing
% pipeline using EEG lab. The data used in this pipeline's development were
% from Jasmine Kwasa's ADHD EEG project, published in "Top-down attention modulates auditory-evoked neural responses in neurotypical,
% but not ADHD, young adults",bioRxiv Feb 2021
% BNR_preprocessing_step1: load in both EEG BDF files for the subject.
% Remove unwanted channels, downsample, and sort into visual and auditory. 
% Output: un-epoched/no ICA just visual and just auditory EEG

% Notes on data format:
% First 120 trials = Visual
% Second 120 trials = Auditory

%%

% Enable EEGLab
addpath('C:\Users\benri\Documents\PhD Year 1\Jasmine EEG Data\eeglab2021.1');

% Add location of bdf files to path
addpath('C:\Users\benri\Documents');

% Ask for dataset name to be used throughout the code
prompt = 'Please enter the name for this dataset, then press Enter: ';
dataset_base_name = input(prompt, 's');

% Start EEGLab
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;


EEG_visual = eeg_emptyset();
EEG_visual.run = [];
EEG_visual.datfile = '';
EEG_auditory = eeg_emptyset();
EEG_auditory.run = [];
EEG_auditory.datfile = '';

for ifile = 1:2 % Loop through both files 

% Import raw BDF file, reference to mastoid channels (65 66), name the file
filename = append('C:\Users\benri\Documents\PhD Year 1\Jasmine EEG Data\BDFs and Results Files\',files_to_load(ifile));
EEG = pop_biosig(char(filename), 'ref',[65 66] ,'refoptions',{'keepref','off'});
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',dataset_base_name,'gui','off');
%EEG = eeg_checkset(EEG); % Check consistency of dataset fields

% Remove unwanted channels (EXG6 EXG7 EXG8), overwrite dataset in memory
EEG = pop_select(EEG, 'nochannel',{'EXG6','EXG7','EXG8'});
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off'); 
%EEG = eeg_checkset(EEG); % Check consistency of dataset fields

% Change the sampling rate, overwrite in memory
EEG = pop_resample( EEG, Fs);
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',dataset_base_name,'overwrite','on','gui','off');

% Use end of trial trigger to divide session into auditory and visual
% NOTE TO SELF: event latencies in the event struct are in samples, and
% times in EEG.times are in ms. 
end_trigger = 3850; % 65290 3850
eventstruct = EEG.event;
end_samples = [];
end_event_indices = [];
all_ievents = [];
for ievent = 1:length(eventstruct)
    if eventstruct(ievent).type == end_trigger
        end_samples(end+1) = eventstruct(ievent).latency;
        end_event_indices(end+1) = eventstruct(ievent).urevent;
        all_ievents(end+1) = ievent;
    end
end
if (length(end_samples) ~= 240)||(length(end_event_indices)~= 240)
    disp('uh oh! Did not find 240 trials! Please manually input the time at which visual trials end:');
    disp(['The time length is: ',num2str(EEG.times(end))])
    pop_eegplot(EEG)
    time_to_split = str2num(input('Enter the time here (in ms):','s')); 
    [~,split_index] = min(abs(EEG.times - time_to_split));

elseif length(end_samples) == 240 && length(end_event_indices) == 240
    split_index = round(end_samples(120)); % index at which to split between visual and auditory
end

%% Put the data and event info into the appropriate structure
if ifile == 1
    EEG_visual1 = pop_select(EEG, 'time',[0 EEG.times(split_index)/1000]);
    EEG_auditory1 = pop_select(EEG, 'time', [EEG.times(split_index)/1000 EEG.times(end)/1000]);
elseif ifile == 2
    EEG_visual2 = pop_select(EEG, 'time',[0 EEG.times(split_index)/1000]);
    EEG_auditory2 = pop_select(EEG, 'time', [EEG.times(split_index)/1000 EEG.times(end)/1000]);
end
% WHAT ELSE DO I NEED TO COPY HERE?

end

% Make new sets for visual and auditory
EEG_visual = pop_mergeset(EEG_visual1, EEG_visual2);
EEG_auditory = pop_mergeset(EEG_auditory1, EEG_auditory2);
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG_visual, CURRENTSET, 'setname',append(dataset_base_name,'Visual Raw'),'gui','off');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG_auditory, CURRENTSET, 'setname',append(dataset_base_name,'Auditory Raw'),'gui','off');

