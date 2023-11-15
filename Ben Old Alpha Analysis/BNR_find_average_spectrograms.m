%% BNR_find_average_spectrograms
%% Define directory
cd /lab_data/barblab/Ben/Oscillations
addpath('/lab_data/barblab/Ben/Oscillations')
addpath('/lab_data/barblab/Ben/Oscillations/Fieldtrip-lite170623')
addpath('/lab_data/barblab/Ben/Oscillations/BNR Alpha Analyses')
addpath('/lab_data/barblab/Ben/Oscillations/BNR Final DATA Files')
addpath('/lab_data/barblab/Ben/Oscillations/eeglab2021.1')

% Script to find the average spectrogram for focal and broad, and run
% statistics
subject_tags = {'018XC','019XX','021XC','022XX','023XX','025MX','026MC','027MC',...
    '028MC','029MX','030MX','031MC',...
    '031XX','032MX','034MX','035MX','039XC','040MC',...
    '051XC','053XX','055XX','061XX','063XX','064XC','067XX',...
    '068XX','070XX','073XX','075XC','076XX','080XX','084MX',...
    '091MC','093XX','096MC','098MX',...
    '100MX','101XC','102XX',...
    '103MX','105MX','110XX',...
    '112MX','200XC'};

%% Find minimum number of correct trials for each condition
for isubject = 1:length(subject_tags)
    curr_subj_tag = split(string(subject_tags(isubject)),'');
    curr_subj_tag = curr_subj_tag(2:6)';
    load(join(string(['AUD_ASA',curr_subj_tag(1:3),'_',curr_subj_tag(4:5),'_1-50_DATA.mat']),''),'SCORE');
    num_correct_trials(isubject) = sum(SCORE.hits(:)==1);
    num_F_N_trials(isubject) = sum((SCORE.F_N(:)'.*SCORE.hits(:)')==1);
    num_F_S1_trials(isubject) = sum((SCORE.F_S1(:)'.*SCORE.hits(:)')==1);
    num_F_S2_trials(isubject) = sum((SCORE.F_S2(:)'.*SCORE.hits(:)')==1);
    num_B_N_trials(isubject) = sum((SCORE.B_N(:)'.*SCORE.hits(:)')==1);
    num_B_S1_trials(isubject) = sum((SCORE.B_S1(:)'.*SCORE.hits(:)')==1);
    num_B_S2_trials(isubject) = sum((SCORE.B_S2(:)'.*SCORE.hits(:)')==1);
    num_focal_trials(isubject) = sum(((SCORE.F_N(:)'.*SCORE.hits(:)')==1) + ((SCORE.F_S1(:)'.*SCORE.hits(:)')==1) + ((SCORE.F_S2(:)'.*SCORE.hits(:)')==1));
    num_broad_trials(isubject) = sum(((SCORE.B_N(:)'.*SCORE.hits(:)')==1) + ((SCORE.B_S1(:)'.*SCORE.hits(:)')==1) + ((SCORE.B_S2(:)'.*SCORE.hits(:)')==1));
end
min_correct_trials = min(num_correct_trials);
min_F_N_trials = min(num_F_N_trials);
min_F_S1_trials = min(num_F_S1_trials);
min_F_S2_trials = min(num_F_S2_trials);
min_B_N_trials = min(num_B_N_trials);
min_B_S1_trials = min(num_B_S1_trials);
min_B_S2_trials = min(num_B_S2_trials);
min_focal_trials = min(num_focal_trials);
min_broad_trials = min(num_broad_trials);

%% For each subject, randomly select G trials from each condition
left_frontocentral_channels = 1:11;
left_parietooccipital_channels = 16:27;
right_frontocentral_channels = 32 + [2:4,7:14];
right_parietooccipital_channels = 32 + [21:32];
channels_to_run = [left_parietooccipital_channels,28:32,right_parietooccipital_channels];
parfor isubject = 1:length(subject_tags)
    curr_subj_tag = split(string(subject_tags(isubject)),'');
    curr_subj_tag = curr_subj_tag(2:6)';
    SCORE_data = load(join(string(['AUD_ASA',curr_subj_tag(1:3),'_',curr_subj_tag(4:5),'_1-50_DATA.mat']),''),'SCORE');
    SCORE = SCORE_data.SCORE;
    POWER_data = load(string(append(subject_tags(isubject),'_spectrogram_bytrial a la Shim.mat')));
    freqs = POWER_data.freqs;
    num_trials_used = POWER_data.num_trials_used;
    power = POWER_data.power;
    times = POWER_data.times;
    curr_size = size(power);
    if curr_size(1) < 240
        power(curr_size(1):240,:,:,:) = nan;
    end
    all_focal_trials = ((SCORE.F_N(:)'.*SCORE.hits(:)')==1) + ((SCORE.F_S1(:)'.*SCORE.hits(:)')==1) + ((SCORE.F_S2(:)'.*SCORE.hits(:)')==1);
    all_focal_trials = find(all_focal_trials == 1);
    all_focal_trials_to_use = all_focal_trials(randi(length(all_focal_trials),1,min_focal_trials));
    all_broad_trials = ((SCORE.B_N(:)'.*SCORE.hits(:)')==1) + ((SCORE.B_S1(:)'.*SCORE.hits(:)')==1) + ((SCORE.B_S2(:)'.*SCORE.hits(:)')==1);
    all_broad_trials = find(all_broad_trials == 1);
    all_broad_trials_to_use = all_broad_trials(randi(length(all_broad_trials),1,min_broad_trials));
    
    all_focal_power(:,:,:,:,isubject) = power(all_focal_trials_to_use,channels_to_run,:,:);
    all_broad_power(:,:,:,:,isubject) = power(all_broad_trials_to_use,channels_to_run,:,:);
end
%%	Average across subject and trial for each of the two conditions
average_focal_power = squeeze(nanmean(all_focal_power,[1,2])); % frequency by time by subject
average_broad_power = squeeze(nanmean(all_broad_power,[1,2])); % frequency by time by subject
save('Average Spectrograms.mat','average_focal_power','average_broad_power')

%% Run permutation test
[clusters, p_values, t_sums, permutation_distribution ] = permutest(average_focal_power,average_broad_power);
save('Spectrogram Permutation Stats.mat','clusters','p_values','t_sums','permutation_distribution')
