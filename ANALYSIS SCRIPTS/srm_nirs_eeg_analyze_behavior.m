%% srm_nirs_eeg_analyze_behavior.m

% Benjamin Richardson
% Created: August 31st, 2023

% Script to analyze behavioral sensitivity (d-prime) for SRM NIRS EEG 1

%BehaviorTable = readtable('/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/SRM-NIRS-EEG Behavior Files/srm-nirs-eeg-1.xlsx','Format','auto');
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
        click_distances_to_remove = find(click_distances < 0.2);
        this_trial_click_times(click_distances_to_remove + 1) = [];

        this_trial_target_all = WordTimesTable(string(WordTimesTable.Condition) == this_trial_masker & WordTimesTable.Run == this_trial_run & string(WordTimesTable.Type) == 'Target',4:end);
        this_trial_target_words = table2array(this_trial_target_all(:,1:2:end));
        this_trial_target_times = table2array(this_trial_target_all(:,2:2:end));

        this_trial_masker_all = WordTimesTable(string(WordTimesTable.Condition) == this_trial_masker & WordTimesTable.Run == this_trial_run & string(WordTimesTable.Type) == 'Masker',4:end);
        this_trial_masker_words = table2array(this_trial_masker_all(:,1:2:end));
        this_trial_masker_times = table2array(this_trial_masker_all(:,2:2:end));

        % Store number of color words in the target and masker
        all_num_target_color_words(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_target_color_words(isubject,string(all_maskers) == string(this_trial_masker)) + sum(ismember(this_trial_target_words,color_words));
        all_num_masker_color_words(isubject,string(all_maskers) == string(this_trial_masker)) = all_num_masker_color_words(isubject,string(all_maskers) == string(this_trial_masker)) + sum(ismember(this_trial_masker_words,color_words));

        % Find just color times in target and masker
        this_trial_target_color_times = this_trial_target_times(ismember(this_trial_target_words,color_words));
        this_trial_target_object_times = this_trial_target_times(~ismember(this_trial_target_words,color_words));
        this_trial_masker_color_times = this_trial_masker_times(ismember(this_trial_masker_words,color_words));

%         disp(["this_trial_condition = ", num2str(this_trial_condition)])
%         disp(["this_trial_masker = ", this_trial_masker])
%         disp(["this trial num masker color words = ", num2str(sum(ismember(this_trial_masker_words,color_words)))])
%         
       %% Hit and False Alarm Windows

       threshold_window_start = 0.4; %0.2
       threshold_window_end =  1.1; % 1.0
       tVec = 1:1/44100:14;
       hit_windows = zeros(1,length(tVec)); % create an empty array to define hit windows
       FA_windows = zeros(1,length(tVec)); % create an empty array to define false alarm windows

        % specify hit windows
        for i = 1:length(this_trial_target_color_times) % for each of the current target color times...
            [~,start_index_hit_window] = min(abs(tVec - (this_trial_target_color_times(i)+threshold_window_start))); % ...the hit window will start threshold_window_start seconds after the word onset
            [~,end_index_hit_window] = min(abs(tVec - (this_trial_target_color_times(i)+threshold_window_end))); % ...the hit window will end threshold_window_end seconds after the word onset

            hit_windows(start_index_hit_window:end_index_hit_window) = 1; % a value of 1 in the vector hit_windows indicate an area where, if a click falls, it will be counted as a hit
        end

        % specify false alarm windows
        for i = 1:length(this_trial_masker_color_times) % for each of the current masker times...
            [~,start_index_FA_window] = min(abs(tVec - (this_trial_masker_color_times(i)+threshold_window_start))); % ...the false alarm window will start threshold_window_start seconds after the word onset
            [~,end_index_FA_window] = min(abs(tVec - (this_trial_masker_color_times(i)+threshold_window_end))); % ...the false alarm window will end threshold_window_end seconds after the word onset

            FA_windows(start_index_FA_window:end_index_FA_window) = 1;
        end

        FA_windows(hit_windows == 1) = 0; % any time there is a hit window, there should not be an FA window

        % ...Calculate the hit rate, FA rate in this trial
        for iclick = 1:length(this_trial_click_times)

            [~,current_click_index] = min(abs(tVec - this_trial_click_times(iclick))); % ...find the time index of that click...

            if hit_windows(current_click_index) == 1 % ...if that click falls within a hit window...
                all_hits(isubject,string(all_maskers) == string(this_trial_masker)) = all_hits(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
            elseif FA_windows(current_click_index) == 1 %...otherwise if that click falls within a false alarm window...
                all_FAs(isubject,string(all_maskers) == string(this_trial_masker)) = all_FAs(isubject,string(all_maskers) == string(this_trial_masker)) + 1;
            else % ...if the click is not counted as either
                clicks_not_counted = clicks_not_counted + 1;
            end
            total_clicks = total_clicks + 1;
        end

        % associate it with the correct condition

        
    end

disp([string(subject_ID(isubject,:)),': ', num2str((clicks_not_counted/total_clicks)*100), '% of clicks not counted'])
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

all_FAs_collapsed_left_and_right = [];
all_FAs_collapsed_left_and_right(1,:) = sum(all_FAs(:,[2,3]),2); % itd50 noise
all_FAs_collapsed_left_and_right(2,:) = sum(all_FAs(:,[10,19]),2); % itd500 noise
all_FAs_collapsed_left_and_right(3,:) = sum(all_FAs(:,[11,17]),2); % ildnat noise
all_FAs_collapsed_left_and_right(4,:) = sum(all_FAs(:,[12,14]),2);% ild10 noise
all_FAs_collapsed_left_and_right(5,:) = sum(all_FAs(:,[5,15]),2); % itd50 speech
all_FAs_collapsed_left_and_right(6,:) = sum(all_FAs(:,[1,7]),2);% itd500 speech
all_FAs_collapsed_left_and_right(7,:) = sum(all_FAs(:,[4,20]),2);% ildnat speech
all_FAs_collapsed_left_and_right(8,:) = sum(all_FAs(:,[6,13]),2);% ild10 speech


all_num_target_color_words_collapsed_left_and_right = [];
all_num_target_color_words_collapsed_left_and_right(1,:) = sum(all_num_target_color_words(:,[2,3]),2);% itd50 noise
all_num_target_color_words_collapsed_left_and_right(2,:) = sum(all_num_target_color_words(:,[10,19]),2); % itd500 noise
all_num_target_color_words_collapsed_left_and_right(3,:) = sum(all_num_target_color_words(:,[11,17]),2);% ildnat noise
all_num_target_color_words_collapsed_left_and_right(4,:) = sum(all_num_target_color_words(:,[12,14]),2);% ild10 noise
all_num_target_color_words_collapsed_left_and_right(5,:) = sum(all_num_target_color_words(:,[5,15]),2);% itd50 speech
all_num_target_color_words_collapsed_left_and_right(6,:) = sum(all_num_target_color_words(:,[1,7]),2);% itd500 speech
all_num_target_color_words_collapsed_left_and_right(7,:) = sum(all_num_target_color_words(:,[4,20]),2);% ildnat speech
all_num_target_color_words_collapsed_left_and_right(8,:) = sum(all_num_target_color_words(:,[6,13]),2);% ild10 speech


all_num_masker_color_words_collapsed_left_and_right = [];
all_num_masker_color_words_collapsed_left_and_right(1,:) = sum(all_num_masker_color_words(:,[2,3]),2);% itd50 noise
all_num_masker_color_words_collapsed_left_and_right(2,:) = sum(all_num_masker_color_words(:,[10,19]),2); % itd500 noise
all_num_masker_color_words_collapsed_left_and_right(3,:) = sum(all_num_masker_color_words(:,[11,17]),2);% ildnat noise
all_num_masker_color_words_collapsed_left_and_right(4,:) = sum(all_num_masker_color_words(:,[12,14]),2);% ild10 noise
all_num_masker_color_words_collapsed_left_and_right(5,:) = sum(all_num_masker_color_words(:,[5,15]),2);% itd50 speech
all_num_masker_color_words_collapsed_left_and_right(6,:) = sum(all_num_masker_color_words(:,[1,7]),2);% itd500 speech
all_num_masker_color_words_collapsed_left_and_right(7,:) = sum(all_num_masker_color_words(:,[4,20]),2);% ildnat speech
all_num_masker_color_words_collapsed_left_and_right(8,:) = sum(all_num_masker_color_words(:,[6,13]),2);% ild10 speech

all_hit_rates = all_hits./all_num_target_color_words;
all_hit_rates_collapsed = all_hits_collapsed_left_and_right./all_num_target_color_words_collapsed_left_and_right;
all_hit_rates_collapsed(all_hit_rates_collapsed == 0) = 0.001;
all_hit_rates_collapsed(all_hit_rates_collapsed == 1) = 0.999;

all_FA_rates = all_FAs./all_num_masker_color_words;
all_FA_rates_collapsed = all_FAs_collapsed_left_and_right./all_num_masker_color_words_collapsed_left_and_right;
all_FA_rates_collapsed(all_FA_rates_collapsed == 0) = 0.001;
%% Hit Rate Figure
figure;boxplot(all_hit_rates_collapsed')
xticks([1:8])
xticklabels({'ITD 50 Noise','ITD 500 Noise','ILD Nat Noise','ILD 10 Noise','ITD 50 Speech','ITD 500 Speech','ILD Nat Speech','ILD 10 Speech'})
ylabel('Hit Rate','FontSize',18)
xlabel('Condition','FontSize',18)

%% D-Prime Figure (Just Speech)
d_primes_speech_masker = norminv(all_hit_rates_collapsed(5:end,:)) - norminv(all_FA_rates_collapsed(5:end,:));
%d_primes_speech_masker(1,:) = d_primes_speech_masker(1,:) - 0.2;
figure;
hold on
for i = 1:size(d_primes_speech_masker,2)
    if d_primes_speech_masker(1,i) > d_primes_speech_masker(2,i) % plot in red
        plot([1:2],d_primes_speech_masker(1:2,i),'Color',[1 0 0])
    else % plot in blue
        plot([1:2],d_primes_speech_masker(1:2,i),'Color',[0 0 1])
    end
    if d_primes_speech_masker(3,i) > d_primes_speech_masker(4,i) % plot in red
       plot([2.5:3.5],d_primes_speech_masker(3:4,i),'Color',[1 0 0])
    else % plot in blue
        plot([2.5:3.5],d_primes_speech_masker(3:4,i),'Color',[0 0 1])
    end
end
e1 = errorbar(1,mean(d_primes_speech_masker(1,:),2,'omitnan'),std(d_primes_speech_masker(1,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e2 = errorbar(2,mean(d_primes_speech_masker(2,:),2,'omitnan'),std(d_primes_speech_masker(2,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e3 = errorbar(2.5,mean(d_primes_speech_masker(3,:),2,'omitnan'),std(d_primes_speech_masker(3,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e4 = errorbar(3.5,mean(d_primes_speech_masker(4,:),2,'omitnan'),std(d_primes_speech_masker(4,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));

e1.Marker = 'o';
e1.MarkerSize = 10;
e1.MarkerFaceColor = 'red';
e1.Color = 'red';
e1.CapSize = 15;
e1.LineWidth = 2;
e2.Marker = 'o';
e2.MarkerSize = 10;
e2.MarkerFaceColor = 'blue';
e2.Color = 'blue';
e2.CapSize = 15;
e2.LineWidth = 2;
e3.Marker = 'o';
e3.MarkerSize = 10;
e3.MarkerFaceColor = 'red';
e3.Color = 'red';
e3.CapSize = 15;
e3.LineWidth = 2;
e4.Marker = 'o';
e4.MarkerSize = 10;
e4.MarkerFaceColor = 'blue';
e4.Color = 'blue';
e4.CapSize = 15;
e4.LineWidth = 2;
ax = gca;
ax.LineWidth = 2;
ax.XTickLabel = cellfun(@(a) ['\bf{' a '}'], ax.XTickLabel, 'UniformOutput',false);

xticks([1,2,2.5,3.5])
xlim([0.75 3.75])
xticklabels({'Small ITD','Large ITD','Natural ILD','Broadband ILD'})
ylabel("Behavioral Sensitivity (d')",'FontSize',18,'FontWeight','bold')
xlabel('Condition','FontSize',18,'FontWeight','bold')
title('Speech Masker Behavior','FontSize',18,'FontWeight','bold')
% % Statistics
% [p,tbl,stats] = anova1(d_primes_speech_masker');
% c = multcompare(stats,'Display','off');
% tbl = array2table(c,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% p_values = c(:,6);
% significant_comparisons = find(p_values < 0.05);
% for i = 1:length(significant_comparisons)
%     xtickvalues = [1,2,2.5,3.5];
%     group_a = c(significant_comparisons(i),1);
%     group_b = c(significant_comparisons(i),2);
%     height = max(d_primes_speech_masker(group_a,:)) + (0.25 - 0.1)*rand(1,1);
%     line([xtickvalues(group_a),xtickvalues(group_b)],[height,height],'Color','k')
%     text(mean([xtickvalues(group_a),xtickvalues(group_b)]),height+0.3,'*','FontSize',24)
% end
%% Hit Rate Figure (Just Noise Masker)
figure;
hold on
for i = 1:size(all_hit_rates_collapsed,2)
    if all_hit_rates_collapsed(1,i) > all_hit_rates_collapsed(2,i) % plot in red
        plot([1:2],all_hit_rates_collapsed(1:2,i),'Color',[1 0 0])
    else % plot in blue
        plot([1:2],all_hit_rates_collapsed(1:2,i),'Color',[0 0 1])
    end
    if all_hit_rates_collapsed(3,i) > all_hit_rates_collapsed(4,i) % plot in red
       plot([2.5:3.5],all_hit_rates_collapsed(3:4,i),'Color',[1 0 0])
    else % plot in blue
        plot([2.5:3.5],all_hit_rates_collapsed(3:4,i),'Color',[0 0 1])
    end
end
e1 = errorbar(1,mean(all_hit_rates_collapsed(1,:),2,'omitnan'),std(all_hit_rates_collapsed(1,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e2 = errorbar(2,mean(all_hit_rates_collapsed(2,:),2,'omitnan'),std(all_hit_rates_collapsed(2,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e3 = errorbar(2.5,mean(all_hit_rates_collapsed(3,:),2,'omitnan'),std(all_hit_rates_collapsed(3,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e4 = errorbar(3.5,mean(all_hit_rates_collapsed(4,:),2,'omitnan'),std(all_hit_rates_collapsed(4,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));

e1.Marker = 'o';
e1.MarkerSize = 10;
e1.MarkerFaceColor = 'red';
e1.Color = 'red';
e1.CapSize = 15;
e1.LineWidth = 2;
e2.Marker = 'o';
e2.MarkerSize = 10;
e2.MarkerFaceColor = 'blue';
e2.Color = 'blue';
e2.CapSize = 15;
e2.LineWidth = 2;
e3.Marker = 'o';
e3.MarkerSize = 10;
e3.MarkerFaceColor = 'red';
e3.Color = 'red';
e3.CapSize = 15;
e3.LineWidth = 2;
e4.Marker = 'o';
e4.MarkerSize = 10;
e4.MarkerFaceColor = 'blue';
e4.Color = 'blue';
e4.CapSize = 15;
e4.LineWidth = 2;
ax = gca;
ax.LineWidth = 2;
ax.XTickLabel = cellfun(@(a) ['\bf{' a '}'], ax.XTickLabel, 'UniformOutput',false);

xticks([1,2,2.5,3.5])
xlim([0.75 3.75])
xticklabels({'Small ITD','Large ITD','Natural ILD','Broadband ILD'})
ylabel('Hit Rate (%)','FontSize',18,'FontWeight','bold')
xlabel('Condition','FontSize',18,'FontWeight','bold')
title('Noise Masker Behavior','FontSize',18)
% % Statistics
% [p,tbl,stats] = anova2(all_hit_rates_collapsed(1:4,:)',1,'off');
% c = multcompare(stats,'Display','off');
% tbl = array2table(c,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% p_values = c(:,6);
% significant_comparisons = find(p_values < 0.05);
% for i = 1:length(significant_comparisons)
%     xtickvalues = [1,2,2.5,3.5];
%     group_a = c(significant_comparisons(i),1);
%     group_b = c(significant_comparisons(i),2);
%     height = max(all_hit_rates_collapsed(group_a,:)) + (0.25 - 0.1)*rand(1,1);
%     line([group_a,group_b],[height,height],'Color','k')
%     text(mean([group_a,group_b]),height+0.3,'*','FontSize',24)
% end

%% Hit Rates (Speech Masker)
figure;
hold on
for i = 1:size(all_hit_rates_collapsed,2)
    if all_hit_rates_collapsed(5,i) > all_hit_rates_collapsed(6,i) % plot in red
        plot([1:2],all_hit_rates_collapsed(5:6,i),'Color',[1 0 0])
    else % plot in blue
        plot([1:2],all_hit_rates_collapsed(5:6,i),'Color',[0 0 1])
    end
    if all_hit_rates_collapsed(7,i) > all_hit_rates_collapsed(8,i) % plot in red
       plot([2.5:3.5],all_hit_rates_collapsed(7:8,i),'Color',[1 0 0])
    else % plot in blue
        plot([2.5:3.5],all_hit_rates_collapsed(7:8,i),'Color',[0 0 1])
    end
end
e1 = errorbar(1,mean(all_hit_rates_collapsed(5,:),2,'omitnan'),std(all_hit_rates_collapsed(5,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e2 = errorbar(2,mean(all_hit_rates_collapsed(6,:),2,'omitnan'),std(all_hit_rates_collapsed(6,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e3 = errorbar(2.5,mean(all_hit_rates_collapsed(7,:),2,'omitnan'),std(all_hit_rates_collapsed(7,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e4 = errorbar(3.5,mean(all_hit_rates_collapsed(8,:),2,'omitnan'),std(all_hit_rates_collapsed(8,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));

e1.Marker = 'o';
e1.MarkerSize = 10;
e1.MarkerFaceColor = 'red';
e1.Color = 'red';
e1.CapSize = 15;
e1.LineWidth = 2;
e2.Marker = 'o';
e2.MarkerSize = 10;
e2.MarkerFaceColor = 'blue';
e2.Color = 'blue';
e2.CapSize = 15;
e2.LineWidth = 2;
e3.Marker = 'o';
e3.MarkerSize = 10;
e3.MarkerFaceColor = 'red';
e3.Color = 'red';
e3.CapSize = 15;
e3.LineWidth = 2;
e4.Marker = 'o';
e4.MarkerSize = 10;
e4.MarkerFaceColor = 'blue';
e4.Color = 'blue';
e4.CapSize = 15;
e4.LineWidth = 2;
ax = gca;
ax.LineWidth = 2;
ax.XTickLabel = cellfun(@(a) ['\bf{' a '}'], ax.XTickLabel, 'UniformOutput',false);

xticks([1,2,2.5,3.5])
xlim([0.75 3.75])
xticklabels({'Small ITD','Large ITD','Natural ILD','Broadband ILD'})
ylabel('Hit Rate (%)','FontSize',18,'FontWeight','bold')
xlabel('Condition','FontSize',18,'FontWeight','bold')
title('Speech Masker Behavior (HIT RATE)','FontSize',18)
% % Statistics
% [p,tbl,stats] = anova2(all_hit_rates_collapsed(1:4,:)',1,'off');
% c = multcompare(stats,'Display','off');
% tbl = array2table(c,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% p_values = c(:,6);
% significant_comparisons = find(p_values < 0.05);
% for i = 1:length(significant_comparisons)
%     xtickvalues = [1,2,2.5,3.5];
%     group_a = c(significant_comparisons(i),1);
%     group_b = c(significant_comparisons(i),2);
%     height = max(all_hit_rates_collapsed(group_a,:)) + (0.25 - 0.1)*rand(1,1);
%     line([group_a,group_b],[height,height],'Color','k')
%     text(mean([group_a,group_b]),height+0.3,'*','FontSize',24)
% end

%% Save data
save('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_Behavior_Results.mat','d_primes_speech_masker','all_hit_rates_collapsed','all_FA_rates_collapsed')
hit_rate_table = array2table(all_hit_rates_collapsed);
writetable(rows2vars(hit_rate_table),'C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_Hit_Rates.csv')
