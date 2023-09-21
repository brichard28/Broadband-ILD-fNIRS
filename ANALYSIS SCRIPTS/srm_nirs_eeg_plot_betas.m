%% srm_nirs_eeg_plot_betas.m

GroupResults = readtable('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\Group Results.csv','Format','auto');
subjects = unique(GroupResults.df);
channels = unique(GroupResults.ch_name);
channels = string(channels);
channels(contains(channels,'hbr')) = [];

conditions = unique(GroupResults.Condition);

% channels
dlpfc_channels = [1,2,3,4,5,6];
stg_channels = 7:14;
left_stg_channels = 11:14;
right_stg_channels = 7:10;

% conditions
attend_left_conditions = 4:2:14;
attend_right_conditions = 5:2:15;
speech_conditions = 10:15;
noise_conditions = 4:9;

all_betas = [];
for isubject = 1:length(subjects)
    for ichannel = 1:length(channels)
        for icondition = 1:length(conditions)
            all_betas(isubject,ichannel,icondition) = GroupResults.theta(GroupResults.df == subjects(isubject) & string(GroupResults.ch_name) == string(channels(ichannel)) & string(GroupResults.Condition) == string(conditions(icondition)) & string(GroupResults.Chroma) == "hbo");
        end
    end
end

% plot speech attend left vs. speech attend right 
all_stg_noise_betas = all_betas(:,stg_channels,4:9);
all_stg_speech_betas = all_betas(:,stg_channels,10:15);



%% INCLUDING ALL CHANNELS
% Speech plot Left Hemisphere (attend left vs. attend right)
upper_ylim = 0.3;
lower_ylim = -0.3;
figure;
subplot(2,2,1)
attend_left_betas = all_betas(:,left_stg_channels,intersect(speech_conditions,attend_left_conditions));
attend_right_betas = all_betas(:,left_stg_channels,intersect(speech_conditions,attend_right_conditions));
x_axis_attend_left = 1:3:7;
x_axis_attend_right = 2:3:87;
p1 = violinplot_T(reshape(attend_left_betas,[size(attend_left_betas,1)*size(attend_left_betas,2),3]),{'1','2','3'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
hold on;p2 =  violinplot_T(reshape(attend_right_betas,[size(attend_right_betas,1)*size(attend_right_betas,2),size(attend_right_betas,3)]),{'1','2','3'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
xticks(1.5:3:16.5)
xticklabels({'ITD 500','ITD 50','ILD 10'})
title('Speech Masker Left STG','FontSize',18)
xlabel('Condition','FontSize',18)
ylabel('Beta','FontSize',18)
ylim([lower_ylim upper_ylim])



% Speech plot Right Hemisphere (attend_left vs. attend_right)
subplot(2,2,2)
attend_left_betas = all_betas(:,right_stg_channels,intersect(speech_conditions,attend_left_conditions));
attend_right_betas = all_betas(:,right_stg_channels,intersect(speech_conditions,attend_right_conditions));
x_axis_attend_left = 1:3:7;
x_axis_attend_right = 2:3:87;
p1 = violinplot_T(reshape(attend_left_betas,[size(attend_left_betas,1)*size(attend_left_betas,2),size(attend_left_betas,3)]),{'1','2','3'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
hold on;p2 =  violinplot_T(reshape(attend_right_betas,[size(attend_right_betas,1)*size(attend_right_betas,2),size(attend_right_betas,3)]),{'1','2','3'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
xticks(1.5:3:16.5)
xticklabels({'ITD 500','ITD 50','ILD 10'})
title('Speech Masker Right STG','FontSize',18)
xlabel('Condition','FontSize',18)
ylabel('Beta','FontSize',18)
ylim([lower_ylim upper_ylim])


% Noise plot Left Hemisphere (attend left vs. attend right)
subplot(2,2,3)
attend_left_betas = all_betas(:,left_stg_channels,intersect(noise_conditions,attend_left_conditions));
attend_right_betas = all_betas(:,left_stg_channels,intersect(noise_conditions,attend_right_conditions));
x_axis_attend_left = 1:3:7;
x_axis_attend_right = 2:3:87;
p1 = violinplot_T(reshape(attend_left_betas,[size(attend_left_betas,1)*size(attend_left_betas,2),size(attend_left_betas,3)]),{'1','2','3'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
hold on;p2 =  violinplot_T(reshape(attend_right_betas,[size(attend_right_betas,1)*size(attend_right_betas,2),size(attend_right_betas,3)]),{'1','2','3'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
xticks(1.5:3:16.5)
xticklabels({'ITD 500','ITD 50','ILD 10'})
title('Noise Masker Left STG','FontSize',18)
xlabel('Condition','FontSize',18)
ylabel('Beta','FontSize',18)
ylim([lower_ylim upper_ylim])


% Noise plot Right hemisphere (attend left vs. attend right)
subplot(2,2,4)
attend_left_betas = all_betas(:,right_stg_channels,intersect(noise_conditions,attend_left_conditions));
attend_right_betas = all_betas(:,right_stg_channels,intersect(noise_conditions,attend_right_conditions));
x_axis_attend_left = 1:3:7;
x_axis_attend_right = 2:3:87;
p1 = violinplot_T(reshape(attend_left_betas,[size(attend_left_betas,1)*size(attend_left_betas,2),size(attend_left_betas,3)]),{'1','2','3'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
hold on;p2 =  violinplot_T(reshape(attend_right_betas,[size(attend_right_betas,1)*size(attend_right_betas,2),size(attend_right_betas,3)]),{'1','2','3'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
xticks(1.5:3:16.5)
xticklabels({'ITD 500','ITD 50','ILD 10'})
title('Noise Masker Right STG','FontSize',18)
xlabel('Condition','FontSize',18)
ylabel('Beta','FontSize',18)
ylim([lower_ylim upper_ylim])


sgtitle('REGULAR DONT TOUCH ALL CHANNELS REPRESENTED HERE')






%% CHOOSING CHANNEL ON EACH SIDE BASED ON ATTEND LEFT or ATTEND RIGHT ITD 500 SPEECH NOISE CONTRAST

itd_500_speech_noise_contrasts_attend_left = all_betas(:,:,10) - all_betas(:,:,4);
itd_500_speech_noise_contrasts_attend_right = all_betas(:,:,11) - all_betas(:,:,5);

itd_500_speech_noise_contrasts_left_stg_attend_left = itd_500_speech_noise_contrasts_attend_left(:,left_stg_channels);
itd_500_speech_noise_contrasts_left_stg_attend_right = itd_500_speech_noise_contrasts_attend_right(:,left_stg_channels);

[~,channels_to_choose_by_subject_left_stg_attend_left] = max(itd_500_speech_noise_contrasts_left_stg_attend_left,[],2);
[~,channels_to_choose_by_subject_left_stg_attend_right] = max(itd_500_speech_noise_contrasts_left_stg_attend_right,[],2);


itd_500_speech_noise_contrasts_right_stg_attend_left = itd_500_speech_noise_contrasts_attend_left(:,right_stg_channels);
itd_500_speech_noise_contrasts_right_stg_attend_right = itd_500_speech_noise_contrasts_attend_right(:,right_stg_channels);

[~,channels_to_choose_by_subject_right_stg_attend_left] = max(itd_500_speech_noise_contrasts_right_stg_attend_left,[],2);
[~,channels_to_choose_by_subject_right_stg_attend_right] = max(itd_500_speech_noise_contrasts_right_stg_attend_right,[],2);

% Speech plot Left Hemisphere (attend left vs. attend right)
figure;
subplot(2,2,1)
attend_left_betas = [];
attend_right_betas = [];
for isubject = 1:length(subjects)
    attend_left_betas(isubject,:,:) = all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg_attend_left(isubject)),intersect(speech_conditions,attend_left_conditions));
    attend_right_betas(isubject,:,:) = all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg_attend_right(isubject)),intersect(speech_conditions,attend_right_conditions));
end
x_axis_attend_left = 1:3:7;
x_axis_attend_right = 2:3:87;
p1 = violinplot_T(attend_left_betas,{'1','2','3'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
hold on;p2 =  violinplot_T(attend_right_betas,{'1','2','3'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
xticks(1.5:3:16.5)
xticklabels({'ITD 500','ITD 50','ILD 10'})
title('Speech Masker Left STG','FontSize',18)
xlabel('Condition','FontSize',18)
ylabel('Beta','FontSize',18)
ylim([lower_ylim upper_ylim])


% Speech plot Right Hemisphere (attend_left vs. attend_right)
subplot(2,2,2)
attend_left_betas = [];
attend_right_betas = [];
for isubject = 1:length(subjects)

attend_left_betas(isubject,:,:) = all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg_attend_left(isubject)),intersect(speech_conditions,attend_left_conditions));
attend_right_betas(isubject,:,:) = all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg_attend_right(isubject)),intersect(speech_conditions,attend_right_conditions));
end
x_axis_attend_left = 1:3:7;
x_axis_attend_right = 2:3:87;
p1 = violinplot_T(attend_left_betas,{'1','2','3'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
hold on;p2 =  violinplot_T(attend_right_betas,{'1','2','3'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
xticks(1.5:3:16.5)
xticklabels({'ITD 500','ITD 50','ILD 10'})
title('Speech Masker Right STG','FontSize',18)
xlabel('Condition','FontSize',18)
ylabel('Beta','FontSize',18)
ylim([lower_ylim upper_ylim])

% Noise plot Left Hemisphere (attend left vs. attend right)
subplot(2,2,3)
attend_left_betas = [];
attend_right_betas = [];
for isubject = 1:length(subjects)

attend_left_betas(isubject,:,:) = all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg_attend_left(isubject)),intersect(noise_conditions,attend_left_conditions));
attend_right_betas(isubject,:,:) = all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg_attend_right(isubject)),intersect(noise_conditions,attend_right_conditions));
end
x_axis_attend_left = 1:3:7;
x_axis_attend_right = 2:3:87;
p1 = violinplot_T(attend_left_betas,{'1','2','3'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
hold on;p2 =  violinplot_T(attend_right_betas,{'1','2','3'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
xticks(1.5:3:16.5)
xticklabels({'ITD 500','ITD 50','ILD 10'})
title('Noise Masker Left STG','FontSize',18)
xlabel('Condition','FontSize',18)
ylabel('Beta','FontSize',18)
ylim([lower_ylim upper_ylim])

% Noise plot Right hemisphere (attend left vs. attend right)
subplot(2,2,4)
attend_left_betas = [];
attend_right_betas = [];
for isubject = 1:length(subjects)

attend_left_betas(isubject,:,:) = all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg_attend_left(isubject)),intersect(noise_conditions,attend_left_conditions));
attend_right_betas(isubject,:,:) = all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg_attend_right(isubject)),intersect(noise_conditions,attend_right_conditions));
end
x_axis_attend_left = 1:3:7;
x_axis_attend_right = 2:3:87;
p1 = violinplot_T(attend_left_betas,{'1','2','3'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
hold on;p2 =  violinplot_T(attend_right_betas,{'1','2','3'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
xticks(1.5:3:16.5)
xticklabels({'ITD 500','ITD 50','ILD 10'})
title('Noise Masker Right STG','FontSize',18)
xlabel('Condition','FontSize',18)
ylabel('Beta','FontSize',18)
ylim([lower_ylim upper_ylim])


sgtitle('Channel Chosen (on each side) for Speech Noise Contrast ITD 500 (Attend Left OR Attend Right, respectively)')