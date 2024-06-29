%% srm_nirs_eeg_1_figure_for_f31

% Figure that shows behavior and block averages for F31 grant

%% Behavior (top two panels)
% load behavior data
user = 'Ben';
if user == 'Ben'
    load('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_Behavior_Results.mat')
elseif user == 'Bon'
    load('/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/SRM-NIRS-EEG-1_Behavior_Results.mat')
end
% ORDER: itd50, itd500, ildnat, ild10

subject_ID = char('NDARVX753BR6','NDARZD647HJ1','NDARBL382XK5','NDARGF569BF3','NDARBA306US5','NDARFD284ZP3','NDARAS648DT4','NDARLM531OY3','NDARXL287BE1','NDARRF358KO3','NDARGT639XS6','NDARFV472HU7','NDARDC882NK4','NDARWB491KR3','NDARNL224RR9','NDARTT639AB1','NDARAZC45TW3','NDARNS784LM2','NDARLB144ZM4','NDARTP382XC8','NDARLJ581GD7','NDARGS283RM9','NDARRED356WS', 'NDARHUG535MO','NDARFIN244AL','NDARKAI888JU','NDARBAA679HA'); % 

num_conditions = 20;

%% Analysis Parameters
method = 'weight'; % 'choose' or 'weight'
mode = 'SPEECH NO BREATH'; % 'SPEECH', 'NOISE', or 'BOTH' (add 'NO BREATH' for no breath)
analysis_type = 'collapsed attend and masker PFC time constant';
area = 'PFC'; % 'PFC' or 'STG'
statistic_to_plot = 'mean'; % 'mean' or 'beta'


%% Plotting Behavior
figure;
hold on
% % ITD Conditions
subplot(2,2,1)
yyaxis left
hold on
plot([1,1.5],[mean(d_primes_speech_masker(1,:),2,'omitnan'),mean(d_primes_speech_masker(2,:),2,'omitnan')],'-b','LineWidth',2)

%plot([1:2],d_primes_speech_masker(1:2,:),'Color',[0.4 0.4 0.4])
e1 = errorbar([1],mean(d_primes_speech_masker(1,:),2,'omitnan'),std(d_primes_speech_masker(1,:),[],2,'omitnan')./(sqrt(size(subject_ID,1)-1)),'Color','b');
e1.Marker = 'o';
e1.MarkerFaceColor = 'w';
e1.MarkerSize = 10;
e1.CapSize = 15;
e1.LineWidth = 2;
e2 = errorbar([1.5],mean(d_primes_speech_masker(2,:),2,'omitnan'),std(d_primes_speech_masker(2,:),[],2,'omitnan')./(sqrt(size(subject_ID,1)-1)),'Color','b');
e2.Marker = 'o';
e2.MarkerFaceColor = 'b';
e2.MarkerSize = 10;
e2.CapSize = 15;
e2.LineWidth = 2;


xlim([0.9 1.6])
ylim([1 3.5])
xticks([1,1.5])


% ILD Conditions
subplot(2,2,3)
yyaxis left
hold on
plot([1,1.5],[mean(d_primes_speech_masker(3,:),2,'omitnan'),mean(d_primes_speech_masker(4,:),2,'omitnan')],'-b','LineWidth',2)

%plot([1:2],d_primes_speech_masker(3:4,:),'Color',[0.4 0.4 0.4])
e3 = errorbar([1],mean(d_primes_speech_masker(3,:),2,'omitnan'),std(d_primes_speech_masker(3,:),[],2,'omitnan')./(sqrt(size(subject_ID,1)-1)),'Color','b');
e3.Marker = 'o';
e3.MarkerFaceColor = 'w';
e3.MarkerSize = 10;
e3.CapSize = 15;
e3.LineWidth = 2;

e4 = errorbar([1.5],mean(d_primes_speech_masker(4,:),2,'omitnan'),std(d_primes_speech_masker(4,:),[],2,'omitnan')./(sqrt(size(subject_ID,1)-1)),'Color','b');
e4.Marker = 'o';
e4.MarkerFaceColor = 'b';
e4.MarkerSize = 10;
e4.CapSize = 15;
e4.LineWidth = 2;


xlim([0.9 1.6])
ylim([1 3.5])
xticks([1,1.5])



%% Block Averages (bottom 2 panels)

if user == 'Ben'
    GroupResults = readtable(append('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\Group Results SRM-NIRS-EEG-1 ',analysis_type,' ',mode, '.csv'),'Format','auto');
    addpath('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\errorbar_files\errorbar_files');
else
    GroupResults = readtable(append('/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/Group Results SRM-NIRS-EEG-1 ',analysis_type,' ',mode, '.csv'),'Format','auto');
    addpath('/home/ben/Documents/GitHub/SRM-NIRS-EEG/errorbar_files/errorbar_files');
end
%subjects = unique(GroupResults.ID);

% channels = unique(GroupResults.ch_name);
% channels = string(channels);
% channels(contains(channels,'hbr')) = [];

conditions = unique(GroupResults.Condition);
conditions(string(conditions) == 'Exhale') = [];
conditions(string(conditions) == 'Inhale') = [];
conditions(string(conditions) == 'Hold') = [];
if contains(analysis_type,'collapsed attend and masker')
    conditions([4 5]) = conditions([5 4]);
end
% channels
dlpfc_channels = [1,2,3,4,5,6]; % 1,2,3,4,5,6
stg_channels = 7:14;
left_stg_channels = 11:14;
right_stg_channels = 7:10;

% time
epoch_time_limits = [-5,25];

all_block_averages = [];
for isubject = 1:size(subject_ID,1)
    
    if user == 'Ben'
        this_subject_table = readtable('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\' + string(subject_ID(isubject,:)) + ' block averages ' + mode +'.csv');
%        this_epochs_deleted = readtable('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\' + string(subject_ID(isubject,:)) + ' epochs ' + mode +'.csv');
    else
        this_subject_table = readtable('/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/' + string(subject_ID(isubject,:)) + ' block averages ' + mode +'.csv');
%        this_epochs_deleted = readtable('/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/' + string(subject_ID(isubject,:)) + ' epochs ' + mode +'.csv');
    end
%     if numel(this_epochs_deleted) > 0
%         this_epochs_deleted = table2array(this_epochs_deleted(2:end,2));
%     else
%         this_epochs_deleted = [];
%     end
    all_epochs = this_subject_table.epoch;
    channels = this_subject_table.Properties.VariableNames(5:end); % ADD IN NANS IF CHANNELS IS LESS THAN IT SHOULD BE
    channels(contains(channels,'Hbr')) = [];
    for ichannel = 1:length(channels)
        for icondition = 1:length(conditions)
            this_condition_epochs = this_subject_table.epoch;
            this_condition_epochs = this_condition_epochs(string(this_subject_table.condition) == conditions(icondition));
            unique_epochs = unique(this_condition_epochs);
            curr_hbo = eval('this_subject_table.' + string(channels(ichannel)));
            z = [];
            for i = 1:length(unique_epochs)
                z(i,:) = curr_hbo(all_epochs == unique_epochs(i));

%                 if ~ismember(i,this_epochs_deleted)% is not rejected
%                     z(i,:) = curr_hbo(all_epochs == unique_epochs(i));
%                 else
%                     z(i,:) = nan*ones(121,1);
%                 end
            end
            all_block_averages(isubject,ichannel,icondition,:) = nanmean(z,1);
        end
    end
end

all_stg_block_averages = all_block_averages(:,stg_channels,:,:);
all_pfc_block_averages = all_block_averages(:,dlpfc_channels,:,:);


% weight each beta value by the inverse of the standard error of the
% GLM fit, and include all channels in the analysis
for isubject = 1:size(subject_ID,1)
    % Store block averages
    stg_block_averages_to_plot(isubject,:,:) = squeeze(nanmean(all_stg_block_averages(isubject,:,:,:),2));
    pfc_block_averages_to_plot(isubject,:,:) = squeeze(nanmean(all_pfc_block_averages(isubject,:,:,:),2));

end

if area == 'PFC'
    block_averages_to_plot = pfc_block_averages_to_plot;
elseif area == 'STG'
    block_averages_to_plot = smoothdata(stg_block_averages_to_plot,3,'SmoothingFactor',0.15);
end

all_betas = [];
all_ses = [];
for isubject = 1:size(subject_ID,1)
    if user == 'Ben'
        this_subject_table = readtable('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\' + string(subject_ID(isubject,:)) + ' block averages ' + mode +'.csv');
 %       this_epochs_deleted = readtable('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\' + string(subject_ID(isubject,:)) + ' epochs ' + mode +'.csv');
    else
        this_subject_table = readtable('/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/' + string(subject_ID(isubject,:)) + ' block averages ' + mode +'.csv');
%        this_epochs_deleted = readtable('/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/' + string(subject_ID(isubject,:)) + ' epochs ' + mode +'.csv');
    end
    channels = this_subject_table.Properties.VariableNames(5:end);
    channels(contains(channels,'Hbr')) = [];
    for ichannel = 1:length(channels)
        this_channel_name = erase(string(channels(ichannel)),'Hbo') + ' hbo';
        for icondition = 1:length(conditions)
            all_betas(isubject,ichannel,icondition) = GroupResults.theta(string(GroupResults.ID) == string(subject_ID(isubject,:)) & string(GroupResults.ch_name) == this_channel_name & string(GroupResults.Condition) == string(conditions(icondition)) & string(GroupResults.Chroma) == "hbo");
            all_ses(isubject,ichannel,icondition) = GroupResults.se(string(GroupResults.ID) == string(subject_ID(isubject,:)) & string(GroupResults.ch_name) == this_channel_name & string(GroupResults.Condition) == string(conditions(icondition)) & string(GroupResults.Chroma) == "hbo");

        end
    end
end

all_stg_betas = all_betas(:,stg_channels,:);
all_pfc_betas = all_betas(:,dlpfc_channels,:);

%% Plot summary statistic on same axes as behavior
time = linspace(epoch_time_limits(1),epoch_time_limits(2),size(all_block_averages,4));
% ITD Conditions
subplot(2,2,1)
hold on
yyaxis right
ylim([0.03 0.13])

[~,time_index_0] = min(abs(time - 4)); %0
[~,time_index_10] = min(abs(time - 7)); %10.8
mean_itd_50 = squeeze(mean(block_averages_to_plot(:,2,time_index_0:time_index_10),3));
mean_itd_500 = squeeze(mean(block_averages_to_plot(:,3,time_index_0:time_index_10),3));
beta_itd_50 = squeeze(mean(all_pfc_betas(:,:,2),2));
beta_itd_500 = squeeze(mean(all_pfc_betas(:,:,3),2));

if statistic_to_plot == 'mean'
    plot([1,1.5],[mean(mean_itd_50),mean(mean_itd_500)],'-r','LineWidth',2)

    e5 = errorbar([1],mean(mean_itd_50),std(mean_itd_50,[],1,'omitnan')./(sqrt(size(subject_ID,1)-1)),'Color','r');
elseif statistic_to_plot == 'beta'
    plot([1,1.5],[mean(beta_itd_50),mean(beta_itd_500)],'-r','LineWidth',2)

    e5 = errorbar([1],mean(beta_itd_50),std(beta_itd_50,[],1,'omitnan')./(sqrt(size(subject_ID,1)-1)),'Color','r');

end
e5.Marker = '^';
e5.MarkerFaceColor = 'w';
e5.MarkerSize = 10;
e5.CapSize = 15;
e5.LineWidth = 2;
if statistic_to_plot == 'mean'
    e6 = errorbar([1.5],mean(mean_itd_500),std(mean_itd_500,[],1,'omitnan')./(sqrt(size(subject_ID,1)-1)),'Color','r');
elseif statistic_to_plot == 'beta'
    e6 = errorbar([1.5],mean(beta_itd_500),std(beta_itd_500,[],1,'omitnan')./(sqrt(size(subject_ID,1)-1)),'Color','r');

end
e6.Marker = '^';
e6.MarkerFaceColor = 'r';
e6.MarkerSize = 10;
e6.CapSize = 15;
e6.LineWidth = 2;


% ILD Conditions
subplot(2,2,3)
hold on
yyaxis right
ylim([0.03 0.14])

mean_ild_70n = squeeze(mean(block_averages_to_plot(:,4,time_index_0:time_index_10),3));
mean_ild_10 = squeeze(mean(block_averages_to_plot(:,5,time_index_0:time_index_10),3));
beta_ild_70n = squeeze(mean(all_pfc_betas(:,:,4),2));
beta_ild_10 = squeeze(mean(all_pfc_betas(:,:,5),2));

if statistic_to_plot == 'mean'
    plot([1,1.5],[mean(mean_ild_70n),mean(mean_ild_10)],'-r','LineWidth',2)

    e7 = errorbar([1],mean(mean_ild_70n),std(mean_ild_70n,[],1,'omitnan')./(sqrt(size(subject_ID,1)-1)),'Color','r','Marker','square');
elseif statistic_to_plot == 'beta'
    plot([1,1.5],[mean(beta_ild_70n),mean(beta_ild_10)],'-r','LineWidth',2)

    e7 = errorbar([1],mean(beta_ild_70n),std(beta_ild_70n,[],1,'omitnan')./(sqrt(size(subject_ID,1)-1)),'Color','r','Marker','square');

end
e7.Marker = 'square';
e7.MarkerFaceColor = 'w';
e7.MarkerSize = 10;
e7.CapSize = 15;
e7.LineWidth = 2;
if statistic_to_plot == 'mean'
e8 = errorbar([1.5],mean(mean_ild_10),std(mean_ild_10,[],1,'omitnan')./(sqrt(size(subject_ID,1)-1)),'Color','r','Marker','square');
elseif statistic_to_plot == 'beta'
e8 = errorbar([1.5],mean(beta_ild_10),std(beta_ild_10,[],1,'omitnan')./(sqrt(size(subject_ID,1)-1)),'Color','r','Marker','square');

end
e8.Marker = 'square';
e8.MarkerFaceColor = 'r';
e8.MarkerSize = 10;
e8.CapSize = 15;
e8.LineWidth = 2;

%% 
all_means = [mean_itd_50,mean_itd_500,mean_ild_70n,mean_ild_10];
all_betas = [beta_itd_50,beta_itd_500,beta_ild_70n,beta_ild_10];
save(['C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\', analysis_type,' ',mode, ' ', area, ' results.mat'],'all_betas','all_means')


%% Plot block averages
%ITD Conditions
ymin = -0.3;
ymax = 0.6;
ymin = -0.125;
ymax = 0.15;
time = linspace(epoch_time_limits(1),epoch_time_limits(2),size(all_block_averages,4));
plotting_fs = 20;
lineprop_list = {'-k',{'or','markerfacecolor',[1,1,1],'MarkerIndices',1:plotting_fs:length(time)},{'or','markerfacecolor',[1,0,0],'MarkerIndices',1:plotting_fs:length(time)},{'or','markerfacecolor',[1,1,1],'MarkerIndices',1:plotting_fs:length(time)},{'or','markerfacecolor',[1,0,0],'MarkerIndices',1:plotting_fs:length(time)}};
% ITD conditions
subplot(2,2,2)
hold on
for icondition = 2:3
    this_lineprop = lineprop_list(icondition);
    shadedErrorBar(time,squeeze(nanmean(block_averages_to_plot(:,icondition,:),1)),squeeze(nanstd(block_averages_to_plot(:,icondition,:),[],1))./(sqrt(size(subject_ID,1)-1)),'lineProps',this_lineprop{1,1});
    hold on
end

legend({'Small ITD','Large ITD'},'FontSize',14,'AutoUpdate','off')
ylim([ymin,ymax])
xline(0,'LineWidth',2)
xline(12.8,'LineWidth',2)

% ILD conditions
subplot(2,2,4)
hold on
for icondition = 4:5
    this_lineprop = lineprop_list(icondition);
    shadedErrorBar(time,squeeze(nanmean(block_averages_to_plot(:,icondition,:),1)),squeeze(nanstd(block_averages_to_plot(:,icondition,:),[],1))./(sqrt(size(subject_ID,1)-1)),'lineProps',this_lineprop{1,1});
    hold on
end
legend({'Natural ILD','Broadband ILD'},'FontSize',14,'AutoUpdate','off')
ylim([ymin,ymax])
xline(0,'LineWidth',2)
xline(12.8,'LineWidth',2)

%% Other axis stuff (figure size, labels, bolding, etc.)
fig = gcf;
fig.Position = [500, 250, 900, 570];

annotation('textbox', [0.01, 1, 0, 0], 'string', 'A','FontSize',36,'FontWeight','bold')
annotation('textbox', [0.45, 1, 0, 0], 'string', 'B','FontSize',36,'FontWeight','bold')

h1 = subplot(2,2,1);
hold on
set(h1, 'Units', 'normalized');
set(h1, 'Position', [0.15, 0.6, 0.2, 0.3]);
set(h1, 'XTickLabelRotationMode','manual');
labelArray = [{'Small','Large'}; {'  ITD ','  ITD '}];
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
xticklabels(tickLabels)
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 14;
ax.FontWeight = 'normal';
yyaxis left
ax = gca;
ax.YColor = 'b';
ylabel({"Behavioral";"Sensitivity (d')"},'FontSize',18,'FontWeight','bold')
yyaxis right
ax = gca;
ax.YColor = 'k';
if statistic_to_plot == 'mean'
    if contains(mode, 'NO BREATH')
        ylabel('Mean \DeltaHbO (\muM)','FontSize',18,'FontWeight','bold')
    else
        ylabel('Mean \DeltaHbO (AU)','FontSize',18,'FontWeight','bold')
    end
elseif statistic_to_plot == 'beta'
    if contains(mode, 'NO BREATH')
        ylabel('Mean Beta (\muM)','FontSize',18,'FontWeight','bold')
    else
        ylabel('Mean Beta (AU)','FontSize',18,'FontWeight','bold')
    end
end
h2 = subplot(2,2,2);
set(h2, 'Units', 'normalized');
set(h2, 'Position', [0.55, 0.6, 0.4, 0.3]);
set(h2, 'XTickLabelRotationMode','manual');
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 14;
ax.FontWeight = 'normal';
xlim(epoch_time_limits)
if contains(mode, 'NO BREATH')
    ylabel('\DeltaHbO (\muM)','FontSize',18,'FontWeight','bold')
else
    ylabel('\DeltaHbO (AU)','FontSize',18,'FontWeight','bold')
end

h3 = subplot(2,2,3);
set(h3, 'Units', 'normalized');
set(h3, 'Position', [0.15, 0.2, 0.2, 0.3]);
set(h3, 'XTickLabelRotationMode','manual');
labelArray = [{'Natural','Broadband'}; {'     ILD ','     ILD '}];
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
xticklabels(tickLabels)
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 14;
ax.FontWeight = 'normal';
xlabel('Condition','FontSize',18,'FontWeight','bold')
yyaxis left
ax = gca;
ax.YColor = 'b';
ylabel({"Behavioral";"Sensitivity (d')"},'FontSize',18,'FontWeight','bold')
yyaxis right
ax = gca;
ax.YColor = 'k';
if statistic_to_plot == 'mean'
    if contains(mode, 'NO BREATH')
        ylabel('Mean \DeltaHbO (\muM)','FontSize',18,'FontWeight','bold')
    else
        ylabel('Mean \DeltaHbO (AU)','FontSize',18,'FontWeight','bold')
    end
elseif statistic_to_plot == 'beta'
    if contains(mode, 'NO BREATH')
        ylabel('Mean Beta (\muM)','FontSize',18,'FontWeight','bold')
    else
        ylabel('Mean Beta (AU)','FontSize',18,'FontWeight','bold')
    end
end

h4 = subplot(2,2,4);
set(h4, 'Units', 'normalized');
set(h4, 'Position', [0.55, 0.2, 0.4, 0.3]);
set(h4, 'XTickLabelRotationMode','manual');
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 14;
ax.FontWeight = 'normal';
xlim(epoch_time_limits)
xlabel('Time re stimulus onset (s)','FontSize',18,'FontWeight','bold')
if contains(mode, 'NO BREATH')
    ylabel('\DeltaHbO (\muM)','FontSize',18,'FontWeight','bold')
else
    ylabel('\DeltaHbO (AU)','FontSize',18,'FontWeight','bold')
end

% [p,tbl,stats] = anova1(all_means);
% multcompare(stats)

%% Separate block averages figure
figure;
lineprop_list = {'-k',{'-r'},{'-b'},{'-r'},{'-b'}};
subplot(1,2,1)
hold on
for icondition = 2:3
    this_lineprop = lineprop_list(icondition);
    plot(time,squeeze(block_averages_to_plot(:,icondition,:)),string(this_lineprop{1,1}));
    hold on
end
title('Small (red) vs. Large (blue) ITD')
subplot(1,2,2)
hold on
for icondition = 4:5
    this_lineprop = lineprop_list(icondition);
    plot(time,squeeze(block_averages_to_plot(:,icondition,:)),string(this_lineprop{1,1}));
    hold on
end
title('Natural (red) vs. Broadband (blue) ILD')