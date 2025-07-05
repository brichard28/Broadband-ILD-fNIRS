%% srm_nirs_eeg_analyze_behavior.m

% Benjamin Richardson
% Created: August 31st, 2023

% Script to analyze behavioral sensitivity (d-prime) for SRM NIRS EEG 1

%BehaviorTable = readtable('/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/SRM-NIRS-EEG Behavior Files/srm-nirs-eeg-1.xlsx','Format','auto');
BehaviorTable = readtable('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG Behavior Files\srm-nirs-eeg-1.xlsx','Format','auto');

subject_ID = char('NDARVX753BR6','NDARZD647HJ1','NDARBL382XK5','NDARGF569BF3','NDARBA306US5',...
                'NDARFD284ZP3','NDARAS648DT4','NDARLM531OY3','NDARXL287BE1','NDARRF358KO3','NDARGT639XS6','NDARDC882NK4',...
                'NDARWB491KR3','NDARNL224RR9','NDARTT639AB1','NDARAZC45TW3',...
                'NDARNS784LM2','NDARLB144ZM4','NDARTP382XC8',...
                'NDARLJ581GD7','NDARGS283RM9','NDARRED356WS','NDARHUG535MO',...
                'NDARFIN244AL','NDARKAI888JU','NDARBAA679HA','NDARUXL573SS',...
                'NDARMOL966PB','NDARGHM426BL','NDARSEW256ZA'); %
num_conditions = 20;
%,
all_hits = zeros(size(subject_ID,1),num_conditions);
all_lead_hits = zeros(size(subject_ID,1),num_conditions);
all_lag_hits = zeros(size(subject_ID,1),num_conditions);

all_FAs = zeros(size(subject_ID,1),num_conditions);
all_lead_FAs = zeros(size(subject_ID,1),num_conditions);
all_lag_FAs = zeros(size(subject_ID,1),num_conditions);

all_objects = zeros(size(subject_ID,1),num_conditions);
all_num_target_color_words = zeros(size(subject_ID,1),num_conditions);
all_num_masker_color_words = zeros(size(subject_ID,1),num_conditions);
all_num_target_object_words = zeros(size(subject_ID,1),num_conditions);

all_num_lead_target_color = zeros(size(subject_ID,1),num_conditions);
all_num_lead_masker_color = zeros(size(subject_ID,1),num_conditions);

all_num_lag_target_color = zeros(size(subject_ID,1),num_conditions);
all_num_lag_masker_color = zeros(size(subject_ID,1),num_conditions);



all_num_hit_windows =  zeros(size(subject_ID,1),num_conditions);
all_num_object_windows =  zeros(size(subject_ID,1),num_conditions);
all_num_FA_windows =  zeros(size(subject_ID,1),num_conditions);

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

rt_fig = figure();

for isubject = 1:size(subject_ID,1) % For each subject...
    disp(subject_ID(isubject,:))
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

    this_subject_color_click_distances = [];
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
        click_distances_to_remove = find(click_distances < 0.2);
        this_trial_click_times(click_distances_to_remove + 1) = [];

        %this_trial_click_times = this_trial_click_times + 0.702;

        this_trial_target_all = WordTimesTable(string(WordTimesTable.Condition) == this_trial_masker & WordTimesTable.Run == this_trial_run & string(WordTimesTable.Type) == 'Target',4:end);
        this_trial_target_words = table2array(this_trial_target_all(:,1:2:end));
        this_trial_target_times = table2array(this_trial_target_all(:,2:2:end));

        this_trial_masker_all = WordTimesTable(string(WordTimesTable.Condition) == this_trial_masker & WordTimesTable.Run == this_trial_run & string(WordTimesTable.Type) == 'Masker',4:end);
        this_trial_masker_words = table2array(this_trial_masker_all(:,1:2:end));
        this_trial_masker_times = table2array(this_trial_masker_all(:,2:2:end));

        this_trial_whether_target_lead = sign(this_trial_masker_times - this_trial_target_times);
        this_trial_whether_target_lead(this_trial_whether_target_lead == -1) = 0; % 0 when masker leads, 1 when target leads


        target_color_lead = this_trial_whether_target_lead(ismember(this_trial_target_words,color_words)); % 0 when color lags (second position within pair), 1 when color leads (first position within pair)
        masker_color_lead = 1 - this_trial_whether_target_lead(ismember(this_trial_masker_words,color_words));


        % Store number of color words in the target and masker
        all_num_target_color_words(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_target_color_words(isubject,string(all_maskers) == string(this_trial_masker)) + sum(ismember(this_trial_target_words,color_words));
        all_num_masker_color_words(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_masker_color_words(isubject,string(all_maskers) == string(this_trial_masker)) + sum(ismember(this_trial_masker_words,color_words));
        all_num_target_object_words(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_target_object_words(isubject,string(all_maskers) == string(this_trial_masker)) + sum(~ismember(this_trial_target_words,color_words));
        
        
        all_num_lead_target_color(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_lead_target_color(isubject,string(all_maskers) == string(this_trial_masker)) + sum(ismember(this_trial_target_words,color_words) & this_trial_whether_target_lead == 1);
        all_num_lag_masker_color(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_lag_masker_color(isubject,string(all_maskers) == string(this_trial_masker)) + sum(ismember(this_trial_masker_words,color_words) & this_trial_whether_target_lead == 1);

        all_num_lag_target_color(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_lag_target_color(isubject,string(all_maskers) == string(this_trial_masker)) + sum(ismember(this_trial_target_words,color_words) & this_trial_whether_target_lead == 0);
        all_num_lead_masker_color(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_lead_masker_color(isubject,string(all_maskers) == string(this_trial_masker)) + sum(ismember(this_trial_masker_words,color_words) & this_trial_whether_target_lead == 0);


        % Find just color times in target and masker
        this_trial_target_color_times = this_trial_target_times(ismember(this_trial_target_words,color_words));
        this_trial_target_object_times = this_trial_target_times(~ismember(this_trial_target_words,color_words));
        this_trial_masker_color_times = this_trial_masker_times(ismember(this_trial_masker_words,color_words));

%         disp(["this_trial_condition = ", num2str(this_trial_condition)])
%         disp(["this_trial_masker = ", this_trial_masker])
%         disp(["this trial num masker color words = ", num2str(sum(ismember(this_trial_masker_words,color_words)))])
%         
       %% Hit and False Alarm Windows

       threshold_window_start = 0.3; %0.2
       threshold_window_end =  1.4; % 1.5
       tVec = 0:1/44100:16;
       hit_windows = zeros(1,length(tVec)); % create an empty array to define hit windows
       lead_hit_windows = zeros(1,length(tVec));
        lag_hit_windows = zeros(1,length(tVec));
       FA_windows = zeros(1,length(tVec)); % create an empty array to define false alarm windows
       object_windows = zeros(1,length(tVec));
       lead_FA_windows = zeros(1,length(tVec));
        lag_FA_windows = zeros(1,length(tVec));

        % specify hit windows
        for i = 1:length(this_trial_target_color_times) % for each of the current target color times...
            [~,start_index_hit_window] = min(abs(tVec - (this_trial_target_color_times(i)+threshold_window_start))); % ...the hit window will start threshold_window_start seconds after the word onset
            [~,end_index_hit_window] = min(abs(tVec - (this_trial_target_color_times(i)+threshold_window_end))); % ...the hit window will end threshold_window_end seconds after the word onset

            hit_windows(start_index_hit_window:end_index_hit_window) = 1; % a value of 1 in the vector hit_windows indicate an area where, if a click falls, it will be counted as a hit
             all_num_hit_windows(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_hit_windows(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
            if target_color_lead(i) == 0 % then target LAGs on this color word
                lag_hit_windows(start_index_hit_window:end_index_hit_window) = 1;
            elseif target_color_lead(i) == 1 % then target LEADS on this trial
                lead_hit_windows(start_index_hit_window:end_index_hit_window) = 1;
            end

        end

        % specify false alarm windows
        for i = 1:length(this_trial_masker_color_times) % for each of the current masker times...
            [~,start_index_FA_window] = min(abs(tVec - (this_trial_masker_color_times(i)+threshold_window_start))); % ...the false alarm window will start threshold_window_start seconds after the word onset
            [~,end_index_FA_window] = min(abs(tVec - (this_trial_masker_color_times(i)+threshold_window_end))); % ...the false alarm window will end threshold_window_end seconds after the word onset

            if any(hit_windows(start_index_FA_window:end_index_FA_window) == 1)
                continue
            else
                FA_windows(start_index_FA_window:end_index_FA_window) = 1;
                all_num_FA_windows(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_FA_windows(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                if masker_color_lead(i) == 0 % then masker LAGs on this trial
                    lag_FA_windows(start_index_FA_window:end_index_FA_window) = 1;
                elseif masker_color_lead(i) == 1 % then masker LEADS on this trial
                    lead_FA_windows(start_index_FA_window:end_index_FA_window) = 1;
                end
            end
        end


        % specify object windows

        for i = 1:length(this_trial_target_object_times)
            [~,start_index_object_window] = min(abs(tVec - (this_trial_target_object_times(i)+threshold_window_start))); % ...the object window will start threshold_window_start seconds after the word onset
            [~,end_index_object_window] = min(abs(tVec - (this_trial_target_object_times(i)+threshold_window_end))); % ...the object window will end threshold_window_end seconds after the word onset

            if any(hit_windows(start_index_object_window:end_index_object_window) == 1)
                continue
            elseif any(FA_windows(start_index_object_window:end_index_object_window) == 1)
                continue
            elseif all(object_windows(start_index_object_window:end_index_object_window) == 0)
                object_windows(start_index_object_window:end_index_object_window) = 1;
                all_num_object_windows(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_object_windows(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
            end
        end

        %         FA_windows(hit_windows == 1) = 0; % any time there is a hit window, there should not be an FA window
        %         object_windows(hit_windows == 1) = 0; % any time there is a hit window, there should not be an object window
        %         object_windows(FA_windows == 1) = 0; % any time there is an FA window, there should not be an object window

        test_vector(itrial,:,:) = FA_windows + hit_windows + object_windows;

        % ...Calculate the hit rate, FA rate in this trial
        for iclick = 1:length(this_trial_click_times)
            [~,current_click_index] = min(abs(tVec - this_trial_click_times(iclick))); % ...find the time index of that click...

            if hit_windows(current_click_index) == 1 % ...if that click falls within a hit window...
                all_hits(isubject,string(all_maskers) == string(this_trial_masker)) = all_hits(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                                % populate lead and lag
                if lead_hit_windows(current_click_index) == 1
                    all_lead_hits(isubject,string(all_maskers) == string(this_trial_masker)) = all_lead_hits(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                elseif lag_hit_windows(current_click_index) == 1
                    all_lag_hits(isubject,string(all_maskers) == string(this_trial_masker)) = all_lag_hits(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                end

            elseif FA_windows(current_click_index) == 1 %...otherwise if that click falls within a false alarm window...
                all_FAs(isubject,string(all_maskers) == string(this_trial_masker)) = all_FAs(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                            % populate lead and lag
                if lead_FA_windows(current_click_index) == 1
                    all_lead_FAs(isubject,string(all_maskers) == string(this_trial_masker)) = all_lead_FAs(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                elseif lag_FA_windows(current_click_index) == 1
                    all_lag_FAs(isubject,string(all_maskers) == string(this_trial_masker)) = all_lag_FAs(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
                end
            elseif object_windows(current_click_index) == 1
                all_objects(isubject,string(all_maskers) == string(this_trial_masker)) = all_objects(isubject,string(all_maskers) == string(this_trial_masker)) + 1;

            else% ...if the click is not counted as either
                clicks_not_counted = clicks_not_counted + 1;
            end
            total_clicks = total_clicks + 1;
        end

        % Calculate the distances to the nearest target color word
        distances_click_to_target_color= [];
        for icolortime = 1:length(this_trial_target_color_times)
            distances_click_to_target_color(icolortime,:) = this_trial_click_times - this_trial_target_color_times(icolortime);
        end
        distances_click_to_target_color(distances_click_to_target_color < 0) = nan;

        %% Find the nearest color time to each click (minimum positive value of click_distances in each column)


        [~,nearest_click] = min(abs(distances_click_to_target_color),[],1); % find the nearest click to each target word
        for i = 1:length(this_trial_click_times)
            if isnan(distances_click_to_target_color(:,i)) == ones(1,length(this_trial_target_color_times)) % all of these clicks were before the first word
                nearest_click(i) = nan;
            else

                this_subject_color_click_distances = [this_subject_color_click_distances, distances_click_to_target_color(nearest_click(i),i)];
            end

        end


        
    end

    save(append(string(subject_ID(isubject,:)),'color_click_distances.mat'),'this_subject_color_click_distances');
    figure(rt_fig)
    subplot(round(size(subject_ID,1)/4),4,isubject)
    histogram(this_subject_color_click_distances,100)
    xlim([0,3])

    [counts,edges] = histcounts(this_subject_color_click_distances,100);
    %[peakValue, peakIndex] = findpeaks(counts); % Find peak value and index
    [peakValue, peakIndex] = max(counts);
    peakXValue = edges(peakIndex); % Get the x-axis value corresponding to the peak



    disp(append(string(subject_ID(isubject,:)),' most frequent color reaction time: ',num2str(peakXValue*1000),' ms '))
    disp(append(string(subject_ID(isubject,:)),': ', num2str((clicks_not_counted/total_clicks)*100), '% of clicks not counted'))

end

%% NEW ORDER = itd50 noise, itd500 noise, ildnat noise, ild10 noise, itd50 speech, itd500 speech, ildnat speech, ild10 speech
all_hits_collapsed_left_and_right = [];
all_hits_collapsed_left_and_right(1,:) = sum(all_hits(:,[2,3]),2); % itd50 noise
all_hits_collapsed_left_and_right(2,:) = sum(all_hits(:,[10,19]),2); % itd500 noise
all_hits_collapsed_left_and_right(3,:) = sum(all_hits(:,[11,17]),2); % ildnat noise
all_hits_collapsed_left_and_right(4,:) = sum(all_hits(:,[12,14]),2); % ild10 noise
all_hits_collapsed_left_and_right(5,:) = sum(all_hits(:,[5,15]),2); % itd50 speech
all_hits_collapsed_left_and_right(6,:) = sum(all_hits(:,[1,7]),2); % itd500 speech
all_hits_collapsed_left_and_right(7,:) = sum(all_hits(:,[4,20]),2); % ildnat speech
all_hits_collapsed_left_and_right(8,:) = sum(all_hits(:,[6,13]),2); % ild10 speech

all_lead_hits_collapsed_left_and_right = [];
all_lead_hits_collapsed_left_and_right(1,:) = sum(all_lead_hits(:,[2,3]),2); % itd50 noise
all_lead_hits_collapsed_left_and_right(2,:) = sum(all_lead_hits(:,[10,19]),2); % itd500 noise
all_lead_hits_collapsed_left_and_right(3,:) = sum(all_lead_hits(:,[11,17]),2); % ildnat noise
all_lead_hits_collapsed_left_and_right(4,:) = sum(all_lead_hits(:,[12,14]),2); % ild10 noise
all_lead_hits_collapsed_left_and_right(5,:) = sum(all_lead_hits(:,[5,15]),2); % itd50 speech
all_lead_hits_collapsed_left_and_right(6,:) = sum(all_lead_hits(:,[1,7]),2); % itd500 speech
all_lead_hits_collapsed_left_and_right(7,:) = sum(all_lead_hits(:,[4,20]),2); % ildnat speech
all_lead_hits_collapsed_left_and_right(8,:) = sum(all_lead_hits(:,[6,13]),2); % ild10 speech

all_lag_hits_collapsed_left_and_right = [];
all_lag_hits_collapsed_left_and_right(1,:) = sum(all_lag_hits(:,[2,3]),2); % itd50 noise
all_lag_hits_collapsed_left_and_right(2,:) = sum(all_lag_hits(:,[10,19]),2); % itd500 noise
all_lag_hits_collapsed_left_and_right(3,:) = sum(all_lag_hits(:,[11,17]),2); % ildnat noise
all_lag_hits_collapsed_left_and_right(4,:) = sum(all_lag_hits(:,[12,14]),2); % ild10 noise
all_lag_hits_collapsed_left_and_right(5,:) = sum(all_lag_hits(:,[5,15]),2); % itd50 speech
all_lag_hits_collapsed_left_and_right(6,:) = sum(all_lag_hits(:,[1,7]),2); % itd500 speech
all_lag_hits_collapsed_left_and_right(7,:) = sum(all_lag_hits(:,[4,20]),2); % ildnat speech
all_lag_hits_collapsed_left_and_right(8,:) = sum(all_lag_hits(:,[6,13]),2); % ild10 speech


all_FAs_collapsed_left_and_right = [];
all_FAs_collapsed_left_and_right(1,:) = sum(all_FAs(:,[2,3]),2); % itd50 noise
all_FAs_collapsed_left_and_right(2,:) = sum(all_FAs(:,[10,19]),2); % itd500 noise
all_FAs_collapsed_left_and_right(3,:) = sum(all_FAs(:,[11,17]),2); % ildnat noise
all_FAs_collapsed_left_and_right(4,:) = sum(all_FAs(:,[12,14]),2);% ild10 noise
all_FAs_collapsed_left_and_right(5,:) = sum(all_FAs(:,[5,15]),2); % itd50 speech
all_FAs_collapsed_left_and_right(6,:) = sum(all_FAs(:,[1,7]),2);% itd500 speech
all_FAs_collapsed_left_and_right(7,:) = sum(all_FAs(:,[4,20]),2);% ildnat speech
all_FAs_collapsed_left_and_right(8,:) = sum(all_FAs(:,[6,13]),2);% ild10 speech

all_lead_FAs_collapsed_left_and_right = [];
all_lead_FAs_collapsed_left_and_right(1,:) = sum(all_lead_FAs(:,[2,3]),2); % itd50 noise
all_lead_FAs_collapsed_left_and_right(2,:) = sum(all_lead_FAs(:,[10,19]),2); % itd500 noise
all_lead_FAs_collapsed_left_and_right(3,:) = sum(all_lead_FAs(:,[11,17]),2); % ildnat noise
all_lead_FAs_collapsed_left_and_right(4,:) = sum(all_lead_FAs(:,[12,14]),2); % ild10 noise
all_lead_FAs_collapsed_left_and_right(5,:) = sum(all_lead_FAs(:,[5,15]),2); % itd50 speech
all_lead_FAs_collapsed_left_and_right(6,:) = sum(all_lead_FAs(:,[1,7]),2); % itd500 speech
all_lead_FAs_collapsed_left_and_right(7,:) = sum(all_lead_FAs(:,[4,20]),2); % ildnat speech
all_lead_FAs_collapsed_left_and_right(8,:) = sum(all_lead_FAs(:,[6,13]),2); % ild10 speech

all_lag_FAs_collapsed_left_and_right = [];
all_lag_FAs_collapsed_left_and_right(1,:) = sum(all_lag_FAs(:,[2,3]),2); % itd50 noise
all_lag_FAs_collapsed_left_and_right(2,:) = sum(all_lag_FAs(:,[10,19]),2); % itd500 noise
all_lag_FAs_collapsed_left_and_right(3,:) = sum(all_lag_FAs(:,[11,17]),2); % ildnat noise
all_lag_FAs_collapsed_left_and_right(4,:) = sum(all_lag_FAs(:,[12,14]),2); % ild10 noise
all_lag_FAs_collapsed_left_and_right(5,:) = sum(all_lag_FAs(:,[5,15]),2); % itd50 speech
all_lag_FAs_collapsed_left_and_right(6,:) = sum(all_lag_FAs(:,[1,7]),2); % itd500 speech
all_lag_FAs_collapsed_left_and_right(7,:) = sum(all_lag_FAs(:,[4,20]),2); % ildnat speech
all_lag_FAs_collapsed_left_and_right(8,:) = sum(all_lag_FAs(:,[6,13]),2); % ild10 speech


all_objects_collapsed_left_and_right = [];
all_objects_collapsed_left_and_right(1,:) = sum(all_objects(:,[2,3]),2); % itd50 noise
all_objects_collapsed_left_and_right(2,:) = sum(all_objects(:,[10,19]),2); % itd500 noise
all_objects_collapsed_left_and_right(3,:) = sum(all_objects(:,[11,17]),2); % ildnat noise
all_objects_collapsed_left_and_right(4,:) = sum(all_objects(:,[12,14]),2);% ild10 noise
all_objects_collapsed_left_and_right(5,:) = sum(all_objects(:,[5,15]),2); % itd50 speech
all_objects_collapsed_left_and_right(6,:) = sum(all_objects(:,[1,7]),2);% itd500 speech
all_objects_collapsed_left_and_right(7,:) = sum(all_objects(:,[4,20]),2);% ildnat speech
all_objects_collapsed_left_and_right(8,:) = sum(all_objects(:,[6,13]),2);% ild10 speech


all_num_target_color_words_collapsed_left_and_right = [];
all_num_target_color_words_collapsed_left_and_right(1,:) = sum(all_num_target_color_words(:,[2,3]),2);% itd50 noise
all_num_target_color_words_collapsed_left_and_right(2,:) = sum(all_num_target_color_words(:,[10,19]),2); % itd500 noise
all_num_target_color_words_collapsed_left_and_right(3,:) = sum(all_num_target_color_words(:,[11,17]),2);% ildnat noise
all_num_target_color_words_collapsed_left_and_right(4,:) = sum(all_num_target_color_words(:,[12,14]),2);% ild10 noise
all_num_target_color_words_collapsed_left_and_right(5,:) = sum(all_num_target_color_words(:,[5,15]),2);% itd50 speech
all_num_target_color_words_collapsed_left_and_right(6,:) = sum(all_num_target_color_words(:,[1,7]),2);% itd500 speech
all_num_target_color_words_collapsed_left_and_right(7,:) = sum(all_num_target_color_words(:,[4,20]),2);% ildnat speech
all_num_target_color_words_collapsed_left_and_right(8,:) = sum(all_num_target_color_words(:,[6,13]),2);% ild10 speech

all_num_lead_target_color_words_collapsed_left_and_right = [];
all_num_lead_target_color_words_collapsed_left_and_right(1,:) = sum(all_num_lead_target_color(:,[2,3]),2);% itd50 noise
all_num_lead_target_color_words_collapsed_left_and_right(2,:) = sum(all_num_lead_target_color(:,[10,19]),2); % itd500 noise
all_num_lead_target_color_words_collapsed_left_and_right(3,:) = sum(all_num_lead_target_color(:,[11,17]),2);% ildnat noise
all_num_lead_target_color_words_collapsed_left_and_right(4,:) = sum(all_num_lead_target_color(:,[12,14]),2);% ild10 noise
all_num_lead_target_color_words_collapsed_left_and_right(5,:) = sum(all_num_lead_target_color(:,[5,15]),2);% itd50 speech
all_num_lead_target_color_words_collapsed_left_and_right(6,:) = sum(all_num_lead_target_color(:,[1,7]),2);% itd500 speech
all_num_lead_target_color_words_collapsed_left_and_right(7,:) = sum(all_num_lead_target_color(:,[4,20]),2);% ildnat speech
all_num_lead_target_color_words_collapsed_left_and_right(8,:) = sum(all_num_lead_target_color(:,[6,13]),2);% ild10 speech


all_num_lag_target_color_words_collapsed_left_and_right = [];
all_num_lag_target_color_words_collapsed_left_and_right(1,:) = sum(all_num_lag_target_color(:,[2,3]),2);% itd50 noise
all_num_lag_target_color_words_collapsed_left_and_right(2,:) = sum(all_num_lag_target_color(:,[10,19]),2); % itd500 noise
all_num_lag_target_color_words_collapsed_left_and_right(3,:) = sum(all_num_lag_target_color(:,[11,17]),2);% ildnat noise
all_num_lag_target_color_words_collapsed_left_and_right(4,:) = sum(all_num_lag_target_color(:,[12,14]),2);% ild10 noise
all_num_lag_target_color_words_collapsed_left_and_right(5,:) = sum(all_num_lag_target_color(:,[5,15]),2);% itd50 speech
all_num_lag_target_color_words_collapsed_left_and_right(6,:) = sum(all_num_lag_target_color(:,[1,7]),2);% itd500 speech
all_num_lag_target_color_words_collapsed_left_and_right(7,:) = sum(all_num_lag_target_color(:,[4,20]),2);% ildnat speech
all_num_lag_target_color_words_collapsed_left_and_right(8,:) = sum(all_num_lag_target_color(:,[6,13]),2);% ild10 speech


all_num_masker_color_words_collapsed_left_and_right = [];
all_num_masker_color_words_collapsed_left_and_right(1,:) = sum(all_num_masker_color_words(:,[2,3]),2);% itd50 noise
all_num_masker_color_words_collapsed_left_and_right(2,:) = sum(all_num_masker_color_words(:,[10,19]),2); % itd500 noise
all_num_masker_color_words_collapsed_left_and_right(3,:) = sum(all_num_masker_color_words(:,[11,17]),2);% ildnat noise
all_num_masker_color_words_collapsed_left_and_right(4,:) = sum(all_num_masker_color_words(:,[12,14]),2);% ild10 noise
all_num_masker_color_words_collapsed_left_and_right(5,:) = sum(all_num_masker_color_words(:,[5,15]),2);% itd50 speech
all_num_masker_color_words_collapsed_left_and_right(6,:) = sum(all_num_masker_color_words(:,[1,7]),2);% itd500 speech
all_num_masker_color_words_collapsed_left_and_right(7,:) = sum(all_num_masker_color_words(:,[4,20]),2);% ildnat speech
all_num_masker_color_words_collapsed_left_and_right(8,:) = sum(all_num_masker_color_words(:,[6,13]),2);% ild10 speech

all_num_lead_masker_color_words_collapsed_left_and_right = [];
all_num_lead_masker_color_words_collapsed_left_and_right(1,:) = sum(all_num_lead_masker_color(:,[2,3]),2);% itd50 noise
all_num_lead_masker_color_words_collapsed_left_and_right(2,:) = sum(all_num_lead_masker_color(:,[10,19]),2); % itd500 noise
all_num_lead_masker_color_words_collapsed_left_and_right(3,:) = sum(all_num_lead_masker_color(:,[11,17]),2);% ildnat noise
all_num_lead_masker_color_words_collapsed_left_and_right(4,:) = sum(all_num_lead_masker_color(:,[12,14]),2);% ild10 noise
all_num_lead_masker_color_words_collapsed_left_and_right(5,:) = sum(all_num_lead_masker_color(:,[5,15]),2);% itd50 speech
all_num_lead_masker_color_words_collapsed_left_and_right(6,:) = sum(all_num_lead_masker_color(:,[1,7]),2);% itd500 speech
all_num_lead_masker_color_words_collapsed_left_and_right(7,:) = sum(all_num_lead_masker_color(:,[4,20]),2);% ildnat speech
all_num_lead_masker_color_words_collapsed_left_and_right(8,:) = sum(all_num_lead_masker_color(:,[6,13]),2);% ild10 speech


all_num_lag_masker_color_words_collapsed_left_and_right = [];
all_num_lag_masker_color_words_collapsed_left_and_right(1,:) = sum(all_num_lag_masker_color(:,[2,3]),2);% itd50 noise
all_num_lag_masker_color_words_collapsed_left_and_right(2,:) = sum(all_num_lag_masker_color(:,[10,19]),2); % itd500 noise
all_num_lag_masker_color_words_collapsed_left_and_right(3,:) = sum(all_num_lag_masker_color(:,[11,17]),2);% ildnat noise
all_num_lag_masker_color_words_collapsed_left_and_right(4,:) = sum(all_num_lag_masker_color(:,[12,14]),2);% ild10 noise
all_num_lag_masker_color_words_collapsed_left_and_right(5,:) = sum(all_num_lag_masker_color(:,[5,15]),2);% itd50 speech
all_num_lag_masker_color_words_collapsed_left_and_right(6,:) = sum(all_num_lag_masker_color(:,[1,7]),2);% itd500 speech
all_num_lag_masker_color_words_collapsed_left_and_right(7,:) = sum(all_num_lag_masker_color(:,[4,20]),2);% ildnat speech
all_num_lag_masker_color_words_collapsed_left_and_right(8,:) = sum(all_num_lag_masker_color(:,[6,13]),2);% ild10 speech


all_num_target_object_words_collapsed_left_and_right = [];
all_num_target_object_words_collapsed_left_and_right(1,:) = sum(all_num_target_object_words(:,[2,3]),2);% itd50 noise
all_num_target_object_words_collapsed_left_and_right(2,:) = sum(all_num_target_object_words(:,[10,19]),2); % itd500 noise
all_num_target_object_words_collapsed_left_and_right(3,:) = sum(all_num_target_object_words(:,[11,17]),2);% ildnat noise
all_num_target_object_words_collapsed_left_and_right(4,:) = sum(all_num_target_object_words(:,[12,14]),2);% ild10 noise
all_num_target_object_words_collapsed_left_and_right(5,:) = sum(all_num_target_object_words(:,[5,15]),2);% itd50 speech
all_num_target_object_words_collapsed_left_and_right(6,:) = sum(all_num_target_object_words(:,[1,7]),2);% itd500 speech
all_num_target_object_words_collapsed_left_and_right(7,:) = sum(all_num_target_object_words(:,[4,20]),2);% ildnat speech
all_num_target_object_words_collapsed_left_and_right(8,:) = sum(all_num_target_object_words(:,[6,13]),2);% ild10 speech

all_num_hit_windows_collapsed_left_and_right = [];
all_num_hit_windows_collapsed_left_and_right(1,:) = sum(all_num_hit_windows(:,[2,3]),2);% itd50 noise
all_num_hit_windows_collapsed_left_and_right(2,:) = sum(all_num_hit_windows(:,[10,19]),2); % itd500 noise
all_num_hit_windows_collapsed_left_and_right(3,:) = sum(all_num_hit_windows(:,[11,17]),2);% ildnat noise
all_num_hit_windows_collapsed_left_and_right(4,:) = sum(all_num_hit_windows(:,[12,14]),2);% ild10 noise
all_num_hit_windows_collapsed_left_and_right(5,:) = sum(all_num_hit_windows(:,[5,15]),2);% itd50 speech
all_num_hit_windows_collapsed_left_and_right(6,:) = sum(all_num_hit_windows(:,[1,7]),2);% itd500 speech
all_num_hit_windows_collapsed_left_and_right(7,:) = sum(all_num_hit_windows(:,[4,20]),2);% ildnat speech
all_num_hit_windows_collapsed_left_and_right(8,:) = sum(all_num_hit_windows(:,[6,13]),2);% ild10 speech

all_num_FA_windows_collapsed_left_and_right = [];
all_num_FA_windows_collapsed_left_and_right(1,:) = sum(all_num_FA_windows(:,[2,3]),2);% itd50 noise
all_num_FA_windows_collapsed_left_and_right(2,:) = sum(all_num_FA_windows(:,[10,19]),2); % itd500 noise
all_num_FA_windows_collapsed_left_and_right(3,:) = sum(all_num_FA_windows(:,[11,17]),2);% ildnat noise
all_num_FA_windows_collapsed_left_and_right(4,:) = sum(all_num_FA_windows(:,[12,14]),2);% ild10 noise
all_num_FA_windows_collapsed_left_and_right(5,:) = sum(all_num_FA_windows(:,[5,15]),2);% itd50 speech
all_num_FA_windows_collapsed_left_and_right(6,:) = sum(all_num_FA_windows(:,[1,7]),2);% itd500 speech
all_num_FA_windows_collapsed_left_and_right(7,:) = sum(all_num_FA_windows(:,[4,20]),2);% ildnat speech
all_num_FA_windows_collapsed_left_and_right(8,:) = sum(all_num_FA_windows(:,[6,13]),2);% ild10 speech

all_num_object_windows_collapsed_left_and_right = [];
all_num_object_windows_collapsed_left_and_right(1,:) = sum(all_num_object_windows(:,[2,3]),2);% itd50 noise
all_num_object_windows_collapsed_left_and_right(2,:) = sum(all_num_object_windows(:,[10,19]),2); % itd500 noise
all_num_object_windows_collapsed_left_and_right(3,:) = sum(all_num_object_windows(:,[11,17]),2);% ildnat noise
all_num_object_windows_collapsed_left_and_right(4,:) = sum(all_num_object_windows(:,[12,14]),2);% ild10 noise
all_num_object_windows_collapsed_left_and_right(5,:) = sum(all_num_object_windows(:,[5,15]),2);% itd50 speech
all_num_object_windows_collapsed_left_and_right(6,:) = sum(all_num_object_windows(:,[1,7]),2);% itd500 speech
all_num_object_windows_collapsed_left_and_right(7,:) = sum(all_num_object_windows(:,[4,20]),2);% ildnat speech
all_num_object_windows_collapsed_left_and_right(8,:) = sum(all_num_object_windows(:,[6,13]),2);% ild10 speech


all_hit_rates = all_hits./all_num_hit_windows;
all_lead_hit_rates = all_lead_hits./all_num_lead_target_color;
all_lag_hit_rates = all_lag_hits./all_num_lag_target_color;

all_hit_rates_collapsed = all_hits_collapsed_left_and_right./all_num_hit_windows_collapsed_left_and_right;
all_lead_hit_rates_collapsed = all_lead_hits_collapsed_left_and_right./all_num_lead_target_color_words_collapsed_left_and_right;
all_lag_hit_rates_collapsed = all_lag_hits_collapsed_left_and_right./all_num_lag_target_color_words_collapsed_left_and_right;


all_hit_rates(all_hit_rates == 0) = 0.001;
all_hit_rates(all_hit_rates >= 1) = 0.999;
all_lead_hit_rates(all_lead_hit_rates == 0) = 0.001;
all_lead_hit_rates(all_lead_hit_rates >= 1) = 0.999;
all_lag_hit_rates(all_lag_hit_rates == 0) = 0.001;
all_lag_hit_rates(all_lag_hit_rates >= 1) = 0.999;

all_hit_rates_collapsed(all_hit_rates_collapsed == 0) = 0.01;
all_hit_rates_collapsed(all_hit_rates_collapsed >= 1) = 0.99;
all_lead_hit_rates_collapsed(all_lead_hit_rates_collapsed == 0) = 0.01;
all_lead_hit_rates_collapsed(all_lead_hit_rates_collapsed >= 1) = 0.99;
all_lag_hit_rates_collapsed(all_lag_hit_rates_collapsed == 0) = 0.01;
all_lag_hit_rates_collapsed(all_lag_hit_rates_collapsed >= 1) = 0.99;

all_FA_rates = all_FAs./all_num_FA_windows;
all_lead_FA_rates = all_lead_FAs./all_num_lead_masker_color;
all_lag_FA_rates = all_lag_FAs./all_num_lag_masker_color;

all_FA_rates_collapsed = all_FAs_collapsed_left_and_right./all_num_FA_windows_collapsed_left_and_right;
all_lead_FA_rates_collapsed = all_lead_FAs_collapsed_left_and_right./all_num_lead_target_color_words_collapsed_left_and_right;
all_lag_FA_rates_collapsed = all_lag_FAs_collapsed_left_and_right./all_num_lag_target_color_words_collapsed_left_and_right;


all_FA_rates(all_FA_rates == 0) = 0.001;
all_FA_rates(all_FA_rates >= 1) = 0.999;
all_lead_FA_rates(all_lead_FA_rates == 0) = 0.001;
all_lead_FA_rates(all_lead_FA_rates >= 1) = 0.999;
all_lag_FA_rates(all_lag_FA_rates == 0) = 0.001;
all_lag_FA_rates(all_lag_FA_rates >= 1) = 0.999;

all_FA_rates_collapsed(all_FA_rates_collapsed == 0) = 0.01;
all_FA_rates_collapsed(all_FA_rates_collapsed >= 1) = 0.99;
all_lead_FA_rates_collapsed(all_lead_FA_rates_collapsed == 0) = 0.01;
all_lead_FA_rates_collapsed(all_lead_FA_rates_collapsed >= 1) = 0.99;
all_lag_FA_rates_collapsed(all_lag_FA_rates_collapsed == 0) = 0.01;
all_lag_FA_rates_collapsed(all_lag_FA_rates_collapsed >= 1) = 0.99;

all_object_rates = all_objects./all_num_object_windows;
all_object_rates_collapsed = all_objects_collapsed_left_and_right./all_num_object_windows_collapsed_left_and_right;
all_object_rates(all_object_rates == 0) = 0.001;
all_object_rates(all_object_rates >= 1) = 0.999;
all_object_rates_collapsed(all_object_rates_collapsed == 0) = 0.01;
all_object_rates_collapsed(all_object_rates_collapsed >= 1) = 0.99;

%% D-prime calculation
d_primes_all = norminv(all_hit_rates) - norminv(all_FA_rates);
d_primes_collapsed = norminv(all_hit_rates_collapsed) - norminv(all_FA_rates_collapsed);
d_primes_speech_masker = norminv(all_hit_rates_collapsed(5:end,:)) - norminv(all_FA_rates_collapsed(5:end,:));

lead_d_primes_collapsed = norminv(all_lead_hit_rates_collapsed) - norminv(all_lead_FA_rates_collapsed);
lag_d_primes_collapsed = norminv(all_lag_hit_rates_collapsed) - norminv(all_lag_FA_rates_collapsed);

%% Save data
save('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_Behavior_Results.mat','d_primes_speech_masker','d_primes_collapsed','all_hit_rates_collapsed','all_FA_rates_collapsed', 'all_object_rates_collapsed','lead_d_primes_collapsed','lag_d_primes_collapsed','all_lead_hit_rates_collapsed','all_lead_FA_rates_collapsed','all_lag_hit_rates_collapsed','all_lag_FA_rates_collapsed')

hit_rate_table = array2table(all_hit_rates_collapsed);
writetable(rows2vars(hit_rate_table),'C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_Hit_Rates.csv')

FA_rate_table = array2table(all_FA_rates_collapsed(5:end,:));
writetable(rows2vars(FA_rate_table),'C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_FA_Rates.csv')

object_rate_table = array2table(all_object_rates_collapsed);
writetable(rows2vars(object_rate_table),'C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_object_Rates.csv')

d_prime_table = array2table(d_primes_collapsed(5:end,:));
writetable(rows2vars(d_prime_table),'C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_d_primes.csv')

lead_hit_rate_table = array2table(all_lead_hit_rates_collapsed);
writetable(rows2vars(lead_hit_rate_table),'C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_Lead_Hit_Rates.csv')

lag_hit_rate_table = array2table(all_lag_hit_rates_collapsed);
writetable(rows2vars(lag_hit_rate_table),'C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_Lag_Hit_Rates.csv')

lead_FA_rate_table = array2table(all_lead_FA_rates_collapsed(5:end,:));
writetable(rows2vars(lead_FA_rate_table),'C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_Lead_FA_Rates.csv')

lag_FA_rate_table = array2table(all_lag_FA_rates_collapsed(5:end,:));
writetable(rows2vars(lag_FA_rate_table),'C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_Lag_FA_Rates.csv')

lead_d_prime_table = array2table(lead_d_primes_collapsed(5:end,:));
writetable(rows2vars(lead_d_prime_table),'C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_Lead_d_primes.csv')

lag_d_prime_table = array2table(lag_d_primes_collapsed(5:end,:));
writetable(rows2vars(lag_d_prime_table),'C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_Lag_d_primes.csv')