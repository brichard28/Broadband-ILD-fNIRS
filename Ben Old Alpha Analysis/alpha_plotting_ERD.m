subject_tags = {'018XC','019XX','021XC','022XX','023XX','025MX','026MC','027MC',...
 '028MC','028XX','029MX','030MX','030XX','031MC',...
'031XX','032MX','032XX','034MX','034XC','035MX','039XC','040MC',...
'051XC','053XX','055XX','061XX','063XX','064XC','067XX',...
'068XX','070XX','073XX','074MC','075XC','076XX','080XX','084MX',...
'084XX','091MC','091XC','093MC','093XX','096MC','096XC','098MX',...
'100MX','101XC','102XX',...
'103MX','103XX','105MX','105XX','110XX',...
'112MX','200XC'};
% 057XC excluded right now because times and freqs are NaN?? But need it in
% there for indexing purposes
%'029XX',


subject_info = readtable('SubjectInfo.xlsx','Sheet','Subject Info');
all_subject_IDs = subject_info.(1); % subject tags from the excel sheet
is_adhd = subject_info.(2); % 1s or 0s for ADHD
is_adhd = is_adhd(ismember(all_subject_IDs,subject_tags));
analysis_info = readtable('SubjectInfo.xlsx','Sheet','Analysis Info');

LFC_channels = 1:11;
LPO_channels = 16:27;
RFC_channels = 32 + [2:4,7:14];
RPO_channels = 32 + [21:32];

channel_select = RPO_channels;% [LPO_channels,28:32,RPO_channels]; %[LFC_channels,33,32+5,32+6,32+15,RFC_channels];
channel_name = 'RPO';

trial_type_select = 1:6; % trial types to plot for ADHD vs NT plot
addpath('BNR Alpha Analyses')
for isubject = 1:length(subject_tags)
    load(string(append(subject_tags(isubject),'_ERD.mat')));
    working_power = power;
    all_powers(isubject,:,:,:) = squeeze(working_power(:,:,channel_select,:)); 
    current_alpha_power = squeeze(mean(all_powers(isubject,:,:,:),3)); % take average over channel
    all_alpha_powers(isubject,:,:) = current_alpha_power;
    [~,timeindex1] = min(abs(times(:)+1600));
    [~,timeindex2] = min(abs(times(:)-0));
    all_alpha_powers_baselined(isubject,:,:) = squeeze(((current_alpha_power - mean(current_alpha_power(:,1:timeindex1),[1,2])))./(std(current_alpha_power(:,1:timeindex1),[],[1,2])));
end

working_alpha = all_alpha_powers_baselined;

trimmed_indices = 3:length(all_alpha_powers_baselined)-2;
% NOTE: need to change to dB
figure;
for isubject = 1:length(subject_tags)
    if is_adhd(isubject) == 1
        p1 = plot(times(trimmed_indices),squeeze(mean(working_alpha(isubject,:,trimmed_indices),2)),'r');
    elseif is_adhd(isubject) == 0
        p2 = plot(times(trimmed_indices),squeeze(mean(working_alpha(isubject,:,trimmed_indices),2)),'k');
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
y = squeeze(mean(working_alpha(is_adhd==1,trial_type_select,trimmed_indices),[1,2]));
y = y';
SEM = std(working_alpha(is_adhd==1,trial_type_select,trimmed_indices),[],[1,2])/sqrt(sum(is_adhd==1)-1);
dy = squeeze(SEM);
dy = dy';
x = times(trimmed_indices);
op = 0.3;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[1 0 0],'linestyle','none','facealpha',op);
p1 = plot(x,squeeze(mean(working_alpha(is_adhd==1,trial_type_select,trimmed_indices),[1,2])),'Color',[1 0 0]);
hold on

y = squeeze(mean(working_alpha(is_adhd==0,trial_type_select,trimmed_indices),[1,2]));
y = y';
SEM = std(working_alpha(is_adhd==0,trial_type_select,trimmed_indices),[],[1,2])/sqrt(sum(is_adhd==0)-1);
dy = squeeze(SEM);
dy = dy';
x = times(trimmed_indices);
op = 0.3;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0 0 1],'linestyle','none','facealpha',op);
p2 = plot(x,squeeze(mean(working_alpha(is_adhd==0,trial_type_select,trimmed_indices),[1,2])),'Color',[0 0 1]);
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
figure;
hold on
hold on
y = squeeze(mean(working_alpha(:,[1,2,3],trimmed_indices),[1,2]));
y= y';
SEM = std(working_alpha(:,[1,2,3],trimmed_indices),[],[1,2])/sqrt(sum(is_adhd)-1); % divide by number of subjects - 1
dy = squeeze(SEM);
dy = dy';
x = times(trimmed_indices);
op = 0.3;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[1 0 0],'linestyle','none','facealpha',op);
p1 = plot(times(trimmed_indices),squeeze(mean(working_alpha(:,[1,2,3],trimmed_indices),[1 2])),'Color',[1 0 0]);
hold on

y = squeeze(mean(working_alpha(:,[4,5,6],trimmed_indices),[1,2]));
y=y';
SEM = std(working_alpha(:,[4,5,6],trimmed_indices),[],[1,2])/sqrt(sum(is_adhd)-1); 
dy = squeeze(SEM);
dy = dy';
x = times(trimmed_indices);
op = 0.3;
fill([x;flipud(x)]',[y-dy;flipud(y+dy)]',[0 0 1],'linestyle','none','facealpha',op);
p2 = plot(times(trimmed_indices),squeeze(mean(working_alpha(:,[4,5,6],trimmed_indices),[1 2])),'Color',[0 0 1]);

%legend({'Focal','Broad'},'AutoUpdate','off');
xlabel('Time (ms)')
ylabel('Power (z-score)')
xline(-1600);
xline(0);
%ylim([-2 1.5])
xlim([times(trimmed_indices(1)),times(trimmed_indices(end))])
title(['All Subjects Focal vs. Broad ',channel_name,' (n = ',num2str(length(is_adhd)),')'])
legend([p1 p2],{'Focal','Broad'})