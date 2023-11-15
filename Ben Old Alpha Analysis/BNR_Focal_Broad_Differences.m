%% BNR_Focal_Broad_Differences

% Does alpha lateralize more in focal or broad?

subject_tags = {'018XC','019XX','021XC','022XX','023XX','025MX','026MC','027MC',...
    '028MC','029MX','030MX','031MC',...
    '031XX','032MX','034MX','035MX','039XC','040MC',...
    '051XC','053XX','055XX','061XX','063XX','064XC','067XX',...
    '068XX','070XX','073XX','075XC','076XX','080XX','084MX',...
    '091MC','093XX','096MC','098MX',...
    '100MX','101XC','102XX',...
    '103MX','105MX','110XX',...
    '112MX','200XC'};

subject_info = readtable('SubjectInfo.xlsx','Sheet','Subject Info');
all_subject_IDs = subject_info.(1); % subject tags from the excel sheet
is_adhd = subject_info.(2); % 1s or 0s for ADHD
is_adhd = is_adhd(ismember(all_subject_IDs,subject_tags));
analysis_info = readtable('SubjectInfo.xlsx','Sheet','Analysis Info');

which_subjects = string(input("Which subjects would you like to test? 'All','ADHD',or'NT': ",'s'));
if which_subjects == 'ADHD'
    subject_tags = subject_tags(is_adhd == 1);
elseif which_subjects == 'NT'
    subject_tags = subject_tags(is_adhd == 0);
elseif which_subjects ~= 'All'
    error('Incorrect input for subject group')
end


LFC_channels = 1:11;
LPO_channels = 16:27;
RFC_channels = 32 + [2:4,7:14];
RPO_channels = 32 + [21:32];

%subject_tags = subject_tags(is_adhd == 1);
focal_or_broad = string(input('focal, broad, or both?: ','s'));
if focal_or_broad == 'focal'
    trial_types = 1:3;
elseif focal_or_broad == 'broad'
    trial_types = 4:6;
elseif focal_or_broad == 'both'
    trial_types = 1:6;
end

addpath('BNR Alpha Analyses')
load('Alpha_Powers_BNR_ 2022_FINAL.mat')
% Get Alpha Power for LPO

channel_select = LPO_channels;
channel_name = 'LPO';

left_alpha = all_alpha_powers_baselined(:,:,channel_select,:);

channel_select = RPO_channels;
channel_name = 'LPO';
right_alpha = all_alpha_powers_baselined(:,:,channel_select,:);

% Calculate max desynchronization on the right and left during the prestim
% period
for isubject = 1:length(subject_tags)
    % max desynch during prestim period
    [~,timeindex1] = min(abs(times + 1600)); % where to start baselining (weird bc of zero padding)
    [~,timeindex2] = min(abs(times + 0));
    [right_alpha_power(isubject),~] = max(abs(mean(right_alpha(isubject,:,:,timeindex1:timeindex2),[2,3])));
    [left_alpha_power(isubject),~] = max(abs(mean(left_alpha(isubject,:,:,timeindex1:timeindex2),[2,3])));

    % sum of power during prestim period
    %     right_alpha_power(isubject) = sum(mean(right_alpha(isubject,:,timeindex1:timeindex2),2));
    %     left_alpha_power(isubject) = sum(mean(left_alpha(isubject,:,timeindex1:timeindex2),2));

    % sum of power after stim onset
    %     right_alpha_power(isubject) = sum(abs(squeeze(mean(right_alpha(isubject,:,timeindex2:end),2))));
    %     left_alpha_power(isubject) = sum(abs(squeeze(mean(left_alpha(isubject,:,timeindex2:end),2))));

    % max of power after stim onset
    % right_alpha_power(isubject) = max(mean(right_alpha(isubject,:,timeindex2:end),2));
    %     left_alpha_power(isubject) = max(mean(left_alpha(isubject,:,timeindex2:end),2));

    % average of power after stim onset
    % right_alpha_power(isubject) = mean(right_alpha(isubject,:,timeindex2:timeindex3),[2,3]);
    % left_alpha_power(isubject) = mean(left_alpha(isubject,:,timeindex2:timeindex3),[2,3]);
    %
end

%figure;histogram(right_alpha_max);hold on;histogram(left_alpha_max)
% figure;histogram(right_alpha_power);hold on;histogram(left_alpha_power)
anova1([right_alpha_power;left_alpha_power]')
xticklabels({'RPO','LPO'})
%anova2([right_alpha_power;left_alpha_power]')

% %% Code for 2 way ANOVA - focal versus broad and left versus right when you have it
% all_left = horzcat(FOCAL_left_alpha_power,BROAD_left_alpha_power);
% all_right = horzcat(FOCAL_right_alpha_power,BROAD_right_alpha_power);
% 
% all_left = all_left(~isoutlier(all_focal));
% all_right = all_right(~isoutlier(all_right));
% 
% figure;histogram(all_left(:),10);hold on;histogram(all_right(:),10)
% legend({'LPO','RPO'})
% xlabel('Maximum of pre-stim alpha "dip"')
% ylabel('Num ocurrences')
% title('LPO vs. RPO maximum of pre-stim alpha "dip"')
% 
% all_focal = horzcat(FOCAL_left_alpha_power,FOCAL_right_alpha_power);
% all_broad = horzcat(BROAD_left_alpha_power,BROAD_left_alpha_power);
% 
% all_focal = all_focal(~isoutlier(all_focal));
% all_broad = all_broad(~isoutlier(all_broad));
% 
% 
% figure;histogram(all_focal(:),10);hold on;histogram(all_broad(:),10)
% legend({'Focal','Broad'})
% xlabel('Maximum of pre-stim alpha "dip"')
% ylabel('Num ocurrences')
% title('Focal vs. Broad maximum of pre-stim alpha "dip"')
% 
% num_reps = length(FOCAL_left_alpha_power); % number of measurements repeated for each case, should equal number of subjects
% % Columns = Left and Right, Rows = Focal and Broad
% data_to_enter = horzcat(all_left',all_right');
% anova2(data_to_enter,num_reps)