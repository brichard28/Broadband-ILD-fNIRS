% subject_tags = {'018XC','019XX','021XC','022XX','023XX','025MX','026MC','027MC',...
%  '028MC','028XX','029MX','030MX','030XX','031MC',...
% '031XX','032MX','032XX','034MX','034XC','035MX','039XC','040MC',...
% '051XC','053XX','055XX','061XX','063XX','064XC','067XX',...
% '068XX','070XX','073XX','074MC','075XC','076XX','080XX','084MX',...
% '084XX','091MC','091XC','093MC','093XX','096MC','096XC','098MX',...
% '100MX','101XC','102XX',...
% '103MX','103XX','105MX','105XX','110XX',...
% '112MX','200XC'};
% 057XC excluded right now because times and freqs are NaN?? But need it in
% there for indexing purposes
%'029XX',


%% Using Only Day 1 Data
subject_tags = {'018XC','019XX','021XC','022XX','023XX','025MX','026MC','027MC',...
    '028MC','029MX','030MX','031MC',...
    '031XX','032MX','034MX','035MX','039XC','040MC',...
    '051XC','053XX','055XX','061XX','063XX','064XC','067XX',...
    '068XX','070XX','073XX','075XC','076XX','080XX','084MX',...
    '091MC','093XX','096MC','098MX',...
    '100MX','101XC','102XX',...
    '103MX','105MX','110XX',...
    '112MX','200XC'};


load('Alpha_Powers_a_la_Shim.mat')
subject_info = readtable('SubjectInfo.xlsx','Sheet','Subject Info');
all_subject_IDs = subject_info.(1); % subject tags from the excel sheet
is_adhd = subject_info.(2); % 1s or 0s for ADHD
is_adhd = is_adhd(ismember(all_subject_IDs,subject_tags));
analysis_info = readtable('SubjectInfo.xlsx','Sheet','Analysis Info');

LFC_channels = 1:11;
LPO_channels = 16:27;
RFC_channels = 32 + [2:4,7:14];
RPO_channels = 32 + [21:32];

channel_select = 1:64;%[LPO_channels,29:32,RPO_channels]; %[LFC_channels,33,32+5,32+6,32+15,RFC_channels];
channel_name = 'APO';

trial_type_select = 4:6; % trial types to plot for ADHD vs NT plot
addpath('C:\Users\benri\Documents\PhD Year 1\Jasmine EEG Data\DATA\BNR Final DATA Files')

working_alpha = all_alpha_powers_baselined;
% for itrial = 1:size(working_alpha,2)
%     if any(isnan(working_alpha(:,itrial,:,:)))
%         working_alpha(:,itrial,:,:) = 0;
%     end
% end

trimmed_indices = 3:size(working_alpha,4)-2;
% NOTE: need to change to dB
figure;
for isubject = 1:length(subject_tags)
    curr_trace = squeeze(mean(working_alpha(isubject,:,channel_select,trimmed_indices),[2,3]));
    if any(isnan(curr_trace)) % remove NaN
        continue
    end
    if is_adhd(isubject) == 1
        p1 = plot(times(trimmed_indices),curr_trace,'r');
    elseif is_adhd(isubject) == 0
        p2 = plot(times(trimmed_indices),curr_trace,'k');
    end
    hold on
end
title(['All individual subjects, all trial types ',channel_name,' (n = ', num2str(length(is_adhd)),')'])
xlabel('Time (ms)')
ylabel('Power (z-score)')
xlim([times(trimmed_indices(1)),times(trimmed_indices(end))])
legend([p1 p2],{'Neurotypical','ADHD'})
%legend({'Neurotypical','ADHD'})

%% Plot ADHD Versus Neurotypical
figure;
hold on
y = squeeze(nanmean(working_alpha(is_adhd==1,:,channel_select,trimmed_indices),[1,2,3]));
SEM = nanstd(working_alpha(is_adhd==1,:,channel_select,trimmed_indices),[],[1,2,3])/sqrt(sum(is_adhd==1)-1);
dy = squeeze(SEM);
x = times(trimmed_indices);
op = 0.3;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[1 0 0],'linestyle','none','FaceAlpha',op);
p1 = plot(x,y,'Color',[1 0 0]);
hold on

y = squeeze(nanmean(working_alpha(is_adhd==0,:,channel_select,trimmed_indices),[1,2,3]));
SEM = nanstd(working_alpha(is_adhd==0,:,channel_select,trimmed_indices),[],[1,2,3])/sqrt(sum(is_adhd==0)-1);
dy = squeeze(SEM);
x = times(trimmed_indices);
op = 0.3;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0 0 1],'linestyle','none','FaceAlpha',op);
p2 = plot(x,y,'Color',[0 0 1]);
%legend({'ADHD','Neurotypical'},'AutoUpdate','off');
xlabel('Time (ms)')
ylabel('Power (z-score)')
xline(-1600);
xline(0);
%ylim([-2 1.5])
xlim([times(trimmed_indices(1)),times(trimmed_indices(end))])
title(['ADHD versus Neurotypical All Trial Types ',channel_name,' (n = ',num2str(sum(is_adhd == 1)),' ADHD, n = ',num2str(sum(is_adhd == 0)),' NT)'])
legend([p1 p2],{'ADHD','Neurotypical'})


%% Plot Focal vs. Broad

focal_alpha_means = [];
broad_alpha_means = [];

[~,time_index_1] = min(abs(times + 0));
[~,time_index_2] = min(abs(times - 1000));

for isubject = 1:length(subject_tags)
    curr_subj_tag = split(string(subject_tags(isubject)),'');
    curr_subj_tag = curr_subj_tag(2:6)';
    load(join(string(['AUD_ASA',curr_subj_tag(1:3),'_',curr_subj_tag(4:5),'_1-50_DATA.mat']),''),'SCORE');
    num_correct_trials(isubject) = sum(SCORE.hits(:)==1);
    focal_trials(isubject,:) = (SCORE.F_N(:)'+SCORE.F_S1(:)' + SCORE.F_S2(:)').*SCORE.hits(:)';
    broad_trials(isubject,:) = (SCORE.B_N(:)'+SCORE.B_S1(:)' + SCORE.B_S2(:)').*SCORE.hits(:)';
    focal_alpha_means(isubject,:) = squeeze(nanmean(working_alpha(isubject,focal_trials(isubject,:)==1,channel_select,time_index_1:time_index_2),[1,2,4]));
    broad_alpha_means(isubject,:) = squeeze(nanmean(working_alpha(isubject,broad_trials(isubject,:)==1,channel_select,time_index_1:time_index_2),[1,2,4]));
    LEFT_focal_alpha_means(isubject,:) = squeeze(nanmean(working_alpha(isubject,focal_trials(isubject,:)==1,LPO_channels,trimmed_indices),[1,2,3]));
    RIGHT_focal_alpha_means(isubject,:) = squeeze(nanmean(working_alpha(isubject,focal_trials(isubject,:)==1,RPO_channels,trimmed_indices),[1,2,3]));
    LEFT_broad_alpha_means(isubject,:) = squeeze(nanmean(working_alpha(isubject,broad_trials(isubject,:)==1,LPO_channels,trimmed_indices),[1,2,3]));
    RIGHT_broad_alpha_means(isubject,:) = squeeze(nanmean(working_alpha(isubject,broad_trials(isubject,:)==1,RPO_channels,trimmed_indices),[1,2,3]));

end

figure;
hold on
y = squeeze(nanmean(focal_alpha_means,1))';
SEM = nanstd(focal_alpha_means,[],1)/sqrt(length(is_adhd)-1); % divide by number of subjects - 1
dy = SEM;
x = times(trimmed_indices);
op = 0.3;
fill([x;flipud(x)],[y-dy';flipud(y+dy')],'k','linestyle','none','FaceAlpha',op);
p1 = plot(x,y,'Color','k');

y = squeeze(nanmean(broad_alpha_means,1))';
SEM = nanstd(broad_alpha_means,[],1)/sqrt(length(is_adhd)-1); % divide by number of subjects - 1
dy = SEM;
x = times(trimmed_indices);
op = 0.3;
fill([x;flipud(x)],[y-dy';flipud(y+dy')],[1 0 0],'linestyle','none','FaceAlpha',op);
p2 = plot(x,y,'Color',[1 0 0]);

%legend({'Focal','Broad'},'AutoUpdate','off');
xlabel('Time (ms)','FontSize',18)
ylabel('Power (z-score)','FontSize',18)
xline(-1600);
xline(0);
xlim([-2500 3200])
%ylim([-2 1.5])
xlim([times(trimmed_indices(1)),times(trimmed_indices(end))])
title(['All Subjects Focal vs. Broad ',channel_name,' (n = ',num2str(length(is_adhd)),')'],'FontSize',18)


%% Stats
[~,stat_index_1] = min(abs(times + 1600));
[~,stat_index_2] = min(abs(times - 3200));
focal_for_stats = focal_alpha_means(:,stat_index_1:stat_index_2);
broad_for_stats = broad_alpha_means(:,stat_index_1:stat_index_2);
time_for_stats = times(stat_index_1:stat_index_2);

[h,p] = ttest(focal_for_stats,broad_for_stats,'Alpha',0.05/length(time_for_stats));
p3 = scatter(time_for_stats(h==1),mean([focal_for_stats(:,h==1);broad_for_stats(:,h==1)]),'c','filled');

legend([p1 p2 p3],{'Focal','Broad','s.d.'})

%% Left and Right on Same plot
figure;
hold on
y = squeeze(nanmean(RIGHT_focal_alpha_means - LEFT_focal_alpha_means,1))';
SEM = nanstd(RIGHT_focal_alpha_means - LEFT_focal_alpha_means,[],1)/sqrt(length(is_adhd)-1); % divide by number of subjects - 1
dy = SEM;
x = times(trimmed_indices);
op = 0.3;
fill([x;flipud(x)],[y-dy';flipud(y+dy')],'m','linestyle','none','FaceAlpha',op);
p1 = plot(x,y,'Color','m');

% y = squeeze(nanmean(RIGHT_focal_alpha_means,1))';
% SEM = nanstd(RIGHT_focal_alpha_means,[],1)/sqrt(length(is_adhd)-1); % divide by number of subjects - 1
% dy = SEM;
% x = times(trimmed_indices);
% op = 0.3;
% fill([x;flipud(x)],[y-dy';flipud(y+dy')],'b','linestyle','none','FaceAlpha',op);
% p3 = plot(x,y,'Color','b');
% 
% xlabel('Time (ms)','FontSize',18)
% ylabel('Power (z-score)','FontSize',18)
% xline(-1600);
% xline(0);
% xlim([-2500 3200])
% %ylim([-2 1.5])
% xlim([times(trimmed_indices(1)),times(trimmed_indices(end))])
% title(['All Subjects FOCAL left vs. right ',channel_name,' (n = ',num2str(length(is_adhd)),')'],'FontSize',18)
% legend([p1 p3],{'Left Focal','Right Focal'})
% 
% 
% figure;
% hold on

y = squeeze(nanmean(RIGHT_broad_alpha_means - LEFT_broad_alpha_means,1))';
SEM = nanstd(RIGHT_broad_alpha_means - LEFT_broad_alpha_means,[],1)/sqrt(length(is_adhd)-1); % divide by number of subjects - 1
dy = SEM;
x = times(trimmed_indices);
op = 0.3;
fill([x;flipud(x)],[y-dy';flipud(y+dy')],[1 0 0],'linestyle','none','FaceAlpha',op);
p2 = plot(x,y,'Color',[1 0 0]);
% 
% y = squeeze(nanmean(RIGHT_broad_alpha_means,1))';
% SEM = nanstd(RIGHT_broad_alpha_means,[],1)/sqrt(length(is_adhd)-1); % divide by number of subjects - 1
% dy = SEM;
% x = times(trimmed_indices);
% op = 0.3;
% fill([x;flipud(x)],[y-dy';flipud(y+dy')],'g','linestyle','none','FaceAlpha',op);
% p4 = plot(x,y,'Color','g');



%legend({'Focal','Broad'},'AutoUpdate','off');
xlabel('Time (ms)','FontSize',18)
ylabel('Power (z-score)','FontSize',18)
xline(-1600);
xline(0);
xlim([-2500 3200])
%ylim([-2 1.5])
xlim([times(trimmed_indices(1)),times(trimmed_indices(end))])
title(['All Subjects Right - Left ',channel_name,' (n = ',num2str(length(is_adhd)),')'],'FontSize',18)

legend([p1 p2],{'Focal','Broad'})

% %% Right minus left plot
% figure;
% hold on
% y = squeeze(nanmean(RIGHT_focal_alpha_means - LEFT_focal_alpha_means,1))';
% SEM = nanstd(RIGHT_focal_alpha_means - LEFT_focal_alpha_means,[],1)/sqrt(length(is_adhd)-1); % divide by number of subjects - 1
% dy = SEM;
% x = times(trimmed_indices);
% op = 0.3;
% fill([x;flipud(x)],[y-dy';flipud(y+dy')],'k','linestyle','none','FaceAlpha',op);
% p1 = plot(x,y,'Color','b');
% y = squeeze(nanmean(RIGHT_broad_alpha_means - LEFT_broad_alpha_means,1))';
% SEM = nanstd(RIGHT_broad_alpha_means - LEFT_broad_alpha_means,[],1)/sqrt(length(is_adhd)-1); % divide by number of subjects - 1
% dy = SEM;
% x = times(trimmed_indices);
% op = 0.3;
% fill([x;flipud(x)],[y-dy';flipud(y+dy')],[1 0 0],'linestyle','none','FaceAlpha',op);
% p2 = plot(x,y,'Color',[1 0 0]);
% xlabel('Time (ms)','FontSize',18)
% ylabel('Power (z-score)','FontSize',18)
% xline(-1600);
% xline(0);
% xlim([-2500 3200])
% %ylim([-2 1.5])
% xlim([times(trimmed_indices(1)),times(trimmed_indices(end))])
% legend([p1 p2],{'Focal Right Minus Left','Broad Right Minus Left'})

% %% Alpha Lateralization Index Plot
% figure;
% hold on
% 
% LEFT_focal_alpha_means = mean(LEFT_focal_alpha_means,1);
% RIGHT_focal_alpha_means = mean(RIGHT_focal_alpha_means,1);
% LEFT_broad_alpha_means = mean(LEFT_broad_alpha_means,1);
% RIGHT_broad_alpha_means = mean(RIGHT_broad_alpha_means,1);
% 
% y = squeeze(nanmean((LEFT_focal_alpha_means - RIGHT_focal_alpha_means)./(LEFT_focal_alpha_means + RIGHT_focal_alpha_means),1))';
% SEM = nanstd((LEFT_focal_alpha_means - RIGHT_focal_alpha_means)./(LEFT_focal_alpha_means + RIGHT_focal_alpha_means),[],1)/sqrt(length(is_adhd)-1); % divide by number of subjects - 1
% y = smoothdata(y,'movmean',100);
% dy = SEM;
% dy = smoothdata(dy,'movmean',100);
% x = times(trimmed_indices);
% op = 0.3;
% fill([x;flipud(x)],[y-dy';flipud(y+dy')],'k','linestyle','none','FaceAlpha',op);
% p1 = plot(x,y,'Color','b');
% y = squeeze(nanmean((LEFT_broad_alpha_means - RIGHT_broad_alpha_means)./(LEFT_broad_alpha_means + RIGHT_broad_alpha_means),1))';
% SEM = nanstd((LEFT_broad_alpha_means - RIGHT_broad_alpha_means)./(LEFT_broad_alpha_means + RIGHT_broad_alpha_means),[],1)/sqrt(length(is_adhd)-1); % divide by number of subjects - 1
% y = smoothdata(y,'movmean',100);
% dy = SEM;
% dy = smoothdata(dy,'movmean',100);
% x = times(trimmed_indices);
% op = 0.3;
% fill([x;flipud(x)],[y-dy';flipud(y+dy')],[1 0 0],'linestyle','none','FaceAlpha',op);
% p2 = plot(x,y,'Color',[1 0 0]);
% xlabel('Time (ms)','FontSize',18)
% ylabel('Alpha Lateralization Index (ALI)','FontSize',18)
% xline(-1600);
% xline(0);
% xlim([-2500 3200])
% %ylim([-2 1.5])
% xlim([times(trimmed_indices(1)),times(trimmed_indices(end))])
% legend([p1 p2],{'Focal ALI','Broad ALI'})