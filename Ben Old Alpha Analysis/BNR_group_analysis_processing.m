%% BNR_group_analysis_processing

% Script to plot group averages of alpha power
% Comparisons:
% For all subjects: All focal versus All Broad
% For neurotypical subjects: All focal versus all broad
% For all ADHD subjects: All focal versus All Broad
% All trial types ADHD vs. Neuro Typical

%% Establish subjects to run and create empty vectors for storage
%'030MX','030XX','031MC',...
%'031XX','032MX','032XX','034MX','034XC','035MX','039XC','040MC',...
%'051XC','053XX','055XX','061XX','063XX','064XC','067XX',...
%'068XX','070XX','073XX','074MC','075XC','076XX','080XX','084MX',...
%'084XX','091MC','091XC','093MC','093XX','096MC','096XC','098MX',...
%'100MX','101XC','102XX','018XC','019XX','021XC','022XX','023XX','025MX','026MC','027MC',...
% '028MC','028XX','029MX','
%'103MX','103XX','105MX','105XX','110XX',...
%'112MX',

% '029XX'
% {'018XC','019XX','021XC','022XX','023XX','025MX','026MC','027MC',...
%     '028MC','028XX','029MX','030MX','030XX','031MC',...
%     '031XX','032MX','032XX','034MX','034XC','035MX','039XC','040MC',...
%     '051XC','053XX','055XX','061XX','063XX','064XC','067XX',...
%     '068XX','070XX','073XX','074MC','075XC','076XX','080XX','084MX',...
%     '084XX','091MC','091XC','093MC','093XX','096MC','096XC','098MX',...
%     '100MX','101XC',
subject_tags = {'018XC','019XX','021XC','022XX','023XX','025MX','026MC','027MC',...
     '028MC','028XX','029MX','030MX','030XX','031MC',...
     '031XX','032MX','032XX','034MX','034XC','035MX','039XC','040MC',...
     '051XC','053XX','055XX','061XX','063XX','064XC','067XX',...
     '068XX','070XX','073XX','074MC','075XC','076XX','080XX','084MX',...
     '084XX','091MC','091XC','093MC','093XX','096MC','096XC','098MX',...
     '100MX','101XC','102XX',...
    '103MX','103XX','105MX','105XX','110XX',...
    '112MX','200XC'};
power_all_subjects = zeros(length(subject_tags),6,64,1792); % subjects x trials x freqs, time

eeglab % need it open for newtimef
fs = 256;
addpath('C:\Users\benri\Documents\PhD Year 1\Jasmine EEG Data\BNR Final DATA Files')
include_incorrect = logical(input('Include incorrect trials? (1 or 0) '));
channel_selection = string(input('Please enter the channels you would like to analyze (form: LPO for Left Parietoocciptal, APO for all parietoocciptal): ','s'));
channel_choices = {'Right FrontoCentral','Left FrontoCentral','Right ParietoOcciptal','Left ParietoOcciptal','All FrontoCentral','All ParietoOcciptal'};
% Left Fronto-Central --> A1-A11
% Left Parieto-Occipital --> A16-A27
% Right Front-Central --> B2-B4,B7-B14,
% Right Parieto-Occipital --> B21-B32

left_frontocentral_channels = 1:11;
left_parietooccipital_channels = 16:27;
right_frontocentral_channels = 32 + [2:4,7:14];
right_parietooccipital_channels = 32 + [21:32];

if channel_selection == 'RFC'
    channels_to_run = right_frontocentral_channels;
    channel_index = 1;
elseif channel_selection == 'LFC'
    channels_to_run = left_frontocentral_channels;
    channel_index = 2;
elseif channel_selection == 'RPO'
    channels_to_run = right_parietooccipital_channels;
    channel_index = 3;
elseif channel_selection == 'LPO'
    channels_to_run = left_parietooccipital_channels;
    channel_index = 4;
elseif channel_selection == 'AFC'
    channels_to_run = [left_frontocentral_channels,right_frontocentral_channels];
    channel_index = 5;
elseif channel_selection == 'APO'
    channels_to_run = [left_parietooccipital_channels,right_parietooccipital_channels];
    channel_index = 6;
elseif channel_selection == 'All'
    channels_to_run = 1:64;
    channel_index = 7;
else
    error('Invalid channel choice, please try again!')
end

%% Other Analysis Choices
subtract_evoked = logical(input('Subtract out mean evoked data? (1 or 0)'));


for isubject = 1:length(subject_tags)
    current_tag = char(subject_tags(isubject));
    num_tag = current_tag(1:3);
    type_tag = current_tag(4:end);
    fstruct = dir(['BNR Final DATA Files\AUD_ASA',num_tag,'_',type_tag,'_1-50_DATA.mat']);
    load(fstruct.name);
    %% Get user input for trial type and channels
    trial_types = {'F_N','F_S1','F_S2','B_N','B_S1','B_S2'};
    type_to_run = 1:6;
    % Specify trials and which trials to run based on input

    all_trials = 1:240;
    if include_incorrect
        F_N_trials = all_trials(logical(SCORE.F_N(:).*~ERP.badTrials')); % F_N_trials
        F_S1_trials = all_trials(logical(SCORE.F_S1(:).*~ERP.badTrials')); % F_S1 trials
        F_S2_trials = all_trials(logical(SCORE.F_S2(:).*~ERP.badTrials')); % F_S2 trials
        B_N_trials = all_trials(logical(SCORE.B_N(:).*~ERP.badTrials')); % B_N trials
        B_S1_trials = all_trials(logical(SCORE.B_S1(:).*~ERP.badTrials')); % B_S1 trials
        B_S2_trials = all_trials(logical(SCORE.B_S2(:).*~ERP.badTrials')); % B_S2 trials
    elseif ~include_incorrect
        F_N_trials = all_trials(logical(SCORE.F_N(:).*SCORE.hits(:).*~ERP.badTrials')); % F_N_trials
        F_S1_trials = all_trials(logical(SCORE.F_S1(:).*SCORE.hits(:).*~ERP.badTrials')); % F_S1 trials
        F_S2_trials = all_trials(logical(SCORE.F_S2(:).*SCORE.hits(:).*~ERP.badTrials')); % F_S2 trials
        B_N_trials = all_trials(logical(SCORE.B_N(:).*SCORE.hits(:).*~ERP.badTrials')); % B_N trials
        B_S1_trials = all_trials(logical(SCORE.B_S1(:).*SCORE.hits(:).*~ERP.badTrials')); % B_S1 trials
        B_S2_trials = all_trials(logical(SCORE.B_S2(:).*SCORE.hits(:).*~ERP.badTrials')); % B_S2 trials
    end

    %% correct for the minimum number of trials for a given condition (we should only use that many trials so we can compare equally)
    all_num_trials = [length(F_N_trials),length(F_S1_trials),length(F_S2_trials),length(B_N_trials),length(B_S1_trials),length(B_S2_trials)];
    min_good_trials = min(all_num_trials);

    %% Further Analysis Parameters
    num_channels = 67;
    Fs = 256;
    epoch_time_limits = [-2.5 4.5];
    timefreq_max_freq = 50; % Hz
    num_wavelet_cycles = 3; % cycles
    padratio = 2; % default from timef function
    timefreq_window_size = (Fs*(epoch_time_limits(2)-epoch_time_limits(1)))/8;
    %freqs = Fs*num_wavelet_cycles/timefreq_window_size*[2:2/padratio:timefreq_window_size]/2; % analysis frequencies as calculated in timef
    %is_inrange = freqs <= timefreq_max_freq; % number of analysis frequencies as calculated in timef (as long as num_wavelet_cycles>0)
    %num_freqs = sum(is_inrange == 1);
    num_times = 200; % default in timef
    for itype = 1:length(type_to_run) % for each trial type...


        %% Preallocate space for calculations
        ersp = []; % num trials X num channels X num frequency bins X num times
        itc = []; % num trials X num channels X num frequency bins X num times
        powbase = []; % should be NaN, we remove the baseline manually
        times = []; % num trials x num channels x num times
        freqs = []; % num trials x num channels x num frequency bins
        erspboot = []; % num trials X num channels X 2 X num times
        itcboot = []; % num trials X num channels X 2 X num times
        num_trials = [];
        ersp_all_trials_all_channels = []; % average over all trials, used to compute mean and std power in the pre-cue period
        ersp_all_trials_only = [];
        itc_all_trials = [];
        powbase_all_trials = [];
        erspboot_all_trials = [];
        itcboot_all_trials = [];
        resting_ersp_mean = [];
        resting_ersp_std = [];
        % Locate the EEG dataset and load it as the variable current_EEG
        current_EEG = ERP.erp;
        % Set the number of channels, number of frames, and the number of trials
        num_frames = size(current_EEG,2);
        num_trials = size(current_EEG,3);
        %% Method 1: Calculation of Spectrogram and grab 8-12 Hz from there
        type = type_to_run(itype);
        trials_to_run = eval(append(string(trial_types(itype)),'_trials'));
        eligible_trials = all_trials(1:num_trials);
        trials_to_run = randsample(eligible_trials,min_good_trials); % Choose the number of trials to run for the selected trial type based on the MINIMUM number of correct/not bad trials for ANY type


        % reset log spectral diffs (ersp), inter-trial phase coherence (itc),
        % baseline power spectrum (powbase), times, freqs, ERSP significance
        % diffs (erspboot), and ITC thresholds (itcboot) for this set
        % For each channel and trial, perform time frequency wavelet analysis
        % using eeglab's timef
        for ichannel = 1:length(channels_to_run)
            for itrial = 1:length(trials_to_run)
                this_EEG = current_EEG(channels_to_run(ichannel),:,trials_to_run(itrial));
                if subtract_evoked % subtract ERP before performing newtimef to get rid of evoked activity
                    this_EEG = this_EEG - mean(current_EEG(channels_to_run(ichannel),:,:),3);
                end
                %% Continuous wavelet transform
                %                 % Zero pad
                %                 n_padsamples = 256*2; % 2 seconds
                %                 eeg_zeropadded = [zeros(1,n_padsamples),this_EEG,zeros(1,n_padsamples)]; % zero padding
                %                 % Calculate the fequency spectrum using wavelet analysis with a
                %                 % Morlet mother wavelet
                %                 adjusted_time_limits = [(epoch_time_limits(1)*1000)-(n_padsamples*1000/fs) (epoch_time_limits(2)*1000)+(n_padsamples*1000/fs)];
                %
                %                 [P,R,PB,t,f,Pboot,Rboot] = newtimef(eeg_zeropadded,length(eeg_zeropadded),adjusted_time_limits,Fs,num_wavelet_cycles,'freqs',[0 timefreq_max_freq],'baseline',NaN,'baseboot',1,'plotersp','off','plotitc','off');%,'subitc','on');
                %
                %                 % We only want the data within the epoch
                %                 [~,timeindex4] = min(abs(t - epoch_time_limits(1)*1000));
                %                 [~,timeindex5] = min(abs(t - epoch_time_limits(2)*1000));
                %                 t = t(timeindex4:timeindex5);
                %                 [~,timeindex3] = min(abs(t + 0));
                %                 [~,timeindex1] = min(abs(t + 2500)); % where to start baselining (weird bc of zero padding)
                %                 [~,timeindex2] = min(abs(t + 1000));
                %                 P = P(:,timeindex4:timeindex5);
                %
                %                 resting_ersp_mean(itrial,ichannel) = mean(P(timeindex1:timeindex2)); % mean during pre-cue period
                %                 resting_ersp_std(itrial,ichannel) = std(P(timeindex1:timeindex2)); % std during pre-cue period
                %
                %                 if resting_ersp_std(itrial,ichannel) == 0
                %                     pause
                %                 end
                %
                %                 ersp(itrial,ichannel,:,:) = P;%(P - resting_ersp_mean(itrial,ichannel))/resting_ersp_std(itrial,ichannel);
                %                 itc(itrial,ichannel,:,:) = R;
                %                 powbase(itrial,ichannel,:,:) = PB;
                %                 times(itrial,ichannel,:) = t;
                %                 freqs(itrial,ichannel,:) = f;
                %                 erspboot(itrial,ichannel,:,:) = Pboot;
                %                 itcboot(itrial,ichannel,:,:) = Rboot;

                %% ERD/ERS
                this_EEG = bandpass(this_EEG,[8 12],256);
                % if normalizing before squaring:
                [~,timeindex1] = min(abs(ERP.t + 0));
                [~,timeindex2] = min(abs(ERP.t + 1600));
                mean_during_precue = mean(this_EEG(1:timeindex2));
                std_during_precue = std(this_EEG(1:timeindex2));
                this_EEG = (this_EEG - mean_during_precue)/std_during_precue;
                this_power(ichannel,:) = (this_EEG).^2;

            end
        end
%         % Average power results over trial and channel
%         ersp_all_trials_all_channels = squeeze(mean(ersp,[1,2]));
%         ersp_all_trials_only = squeeze(mean(ersp,1));
%         itc_all_trials = squeeze(mean(itc,[1,2]));
%         powbase_all_trials = squeeze(mean(powbase,[1,2]));
%         erspboot_all_trials = squeeze(mean(erspboot,[1,2]));
%         itcboot_all_trials = squeeze(mean(itcboot,[1,2]));
% 
%         times = squeeze(mean(times,[1,2]));
%         freqs = squeeze(mean(freqs,[1,2]));

        % Save into relevant vectors
        power_all_subjects(isubject,itype,:,:,:) = this_power; % subjects x trial types X channels X freqs X times
    end
    power = power_all_subjects(isubject,:,:,:);
    filename = append(string(subject_tags(isubject)),'_ERD.mat');
    num_trials_used = length(trials_to_run);
    times = ERP.t;
    save(filename,'power','num_trials_used','times');

end

pause