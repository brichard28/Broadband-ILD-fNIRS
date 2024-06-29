BehaviorTable = readtable('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG Behavior Files\srm-nirs-eeg-1.xlsx','Format','auto');

subject_ID = char('NDARVX753BR6', 'NDARZD647HJ1', 'NDARBL382XK5', 'NDARGF569BF3', 'NDARBA306US5', 'NDARFD284ZP3',...
                   'NDARAS648DT4','NDARLM531OY3', 'NDARXL287BE1', 'NDARRF358KO3', 'NDARGT639XS6', 'NDARDC882NK4',...
                   'NDARWB491KR3','NDARNL224RR9', 'NDARTT639AB1', 'NDARAZC45TW3', 'NDARNS784LM2', 'NDARLB144ZM4', 'NDARTP382XC8',...
                   'NDARLJ581GD7','NDARGS283RM9', 'NDARRED356WS', 'NDARHUG535MO','NDARFIN244AL',...
                   'NDARKAI888JU','NDARBAA679HA'); %
num_conditions = 20;
%,
all_hits = zeros(size(subject_ID,1),num_conditions);
all_FAs = zeros(size(subject_ID,1),num_conditions);
all_num_target_color_words = zeros(size(subject_ID,1),num_conditions);
all_num_masker_color_words = zeros(size(subject_ID,1),num_conditions);

all_maskers = {'m_speech__ild_0__itd_500__targ_r__control_0',...
'm_noise__ild_0__itd_50__targ_l__control_0',...
'm_noise__ild_0__itd_50__targ_r__control_0',...
'm_speech__ild_70n__itd_0__targ_r__control_0',...
'm_speech__ild_0__itd_50__targ_l__control_0',...
'm_speech__ild_10__itd_0__targ_l__control_0',...
'm_speech__ild_0__itd_500__targ_l__control_0',...
'm_speech__ild_0__itd_500__targ_l__control_1',...
'm_noise__ild_0__itd_500__targ_r__control_1',...
'm_noise__ild_0__itd_500__targ_r__control_0',...
'm_noise__ild_70n__itd_0__targ_l__control_0',...
'm_noise__ild_10__itd_0__targ_l__control_0',...
'm_speech__ild_10__itd_0__targ_r__control_0',...
'm_noise__ild_10__itd_0__targ_r__control_0',...
'm_speech__ild_0__itd_50__targ_r__control_0',...
'm_speech__ild_0__itd_500__targ_r__control_1',...
'm_noise__ild_70n__itd_0__targ_r__control_0',...
'm_noise__ild_0__itd_500__targ_l__control_1',...
'm_noise__ild_0__itd_500__targ_l__control_0',...
'm_speech__ild_70n__itd_0__targ_l__control_0'}; % we will maintain this order throughout


all_target_distances = [];
all_masker_distances = [];

for isubject = 1:size(subject_ID,1) % For each subject...
    clicks_not_counted = 0;
    total_clicks = 0;

    % Load the word times for this subject
    %WordTimesTable = readtable("/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/SRM-NIRS-EEG Behavior Files/srm-nirs-eeg-1__s_" + string(subject_ID(isubject,:)) + "__Word_Times.csv");
    WordTimesTable = readtable("C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG Behavior Files\srm-nirs-eeg-1__s_" + string(subject_ID(isubject,:)) + "__Word_Times.csv");

    run_count_per_condition = -1*ones(1,num_conditions); % array to keep track of which run in each condition we are on

    % Find the rows associated with this subject
    rows_this_subject = find(BehaviorTable.S == string(subject_ID(isubject,:)));

    % define color and object words
    color_words = ["red","white","green","blue"];
    object_words = [];

    % For each trial....
    for itrial = 1:length(rows_this_subject)

        this_trial_condition = BehaviorTable.Condition(rows_this_subject(itrial)); % find the condition for this trial
        this_trial_masker = BehaviorTable.masker(rows_this_subject(itrial)); % find the masker type for this trial
        run_count_per_condition(string(all_maskers) == string(this_trial_masker)) = run_count_per_condition(string(all_maskers) == string(this_trial_masker)) + 1;

        this_trial_run = run_count_per_condition(string(all_maskers) == string(this_trial_masker)); % find how many runs of this condition have happened already
        this_trial_click_times = table2array(BehaviorTable(rows_this_subject(itrial),9:end)); % find the click times for this trial
        this_trial_click_times(isnan(this_trial_click_times)) = []; % remove NaN from these click times
        % remove double clicks
        click_distances = diff(this_trial_click_times);
        click_distances_to_remove = find(click_distances < 0.5);
        this_trial_click_times(click_distances_to_remove + 1) = [];

        this_trial_target_all = WordTimesTable(string(WordTimesTable.Condition) == this_trial_masker & WordTimesTable.Run == this_trial_run & string(WordTimesTable.Type) == 'Target',4:end);
        this_trial_target_words = table2array(this_trial_target_all(:,1:2:end));
        this_trial_target_times = table2array(this_trial_target_all(:,2:2:end));

        this_trial_masker_all = WordTimesTable(string(WordTimesTable.Condition) == this_trial_masker & WordTimesTable.Run == this_trial_run & string(WordTimesTable.Type) == 'Masker',4:end);
        this_trial_masker_words = table2array(this_trial_masker_all(:,1:2:end));
        this_trial_masker_times = table2array(this_trial_masker_all(:,2:2:end));

                % Find just color times in target and masker
        this_trial_target_color_times = this_trial_target_times(ismember(this_trial_target_words,color_words));
        this_trial_target_object_times = this_trial_target_times(~ismember(this_trial_target_words,color_words));
        this_trial_masker_color_times = this_trial_masker_times(ismember(this_trial_masker_words,color_words));

        for iclick = 1:length(this_trial_click_times)
            if ~isempty(this_trial_target_color_times)
                distance_to_targets = this_trial_click_times(iclick) - this_trial_target_color_times ;
                distance_to_targets(distance_to_targets < 0) = [];
                all_target_distances = [all_target_distances, distance_to_targets];
            end

            if ~isempty(this_trial_masker_color_times)
                distance_to_maskers = this_trial_click_times(iclick) - this_trial_masker_color_times;
                distance_to_maskers(distance_to_maskers < 0) = [];

                all_masker_distances = [all_masker_distances, distance_to_maskers];
            end

        end
    end
end

figure;
histogram(all_target_distances,400);
xticks(0:0.25:14)
xlabel('Time (s)')
title('Distance from click to nearest preceding target color word','FontSize',18)
ylabel('Frequency of occurrence (count)','FontSize',18)

figure;
histogram(all_masker_distances,400);
xticks(0:0.25:14)
xlabel('Time (s)')
title('Distance from click to nearest preceding masker color word','FontSize',18)
ylabel('Frequency of occurrence (count)','FontSize',18)