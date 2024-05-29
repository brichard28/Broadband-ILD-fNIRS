%% srm_nirs_eeg_plot_speech_vs_noise.m

% version of srm_nirs_eeg_plot_nirs_results.m for just speech vs noise

user = 'Ben';

subject_ID = char('NDARVX753BR6','NDARZD647HJ1','NDARBL382XK5','NDARGF569BF3','NDARBA306US5','NDARFD284ZP3','NDARAS648DT4','NDARLM531OY3','NDARXL287BE1','NDARRF358KO3','NDARGT639XS6','NDARDC882NK4','NDARWB491KR3','NDARNL224RR9','NDARTT639AB1','NDARAZC45TW3','NDARNS784LM2','NDARLB144ZM4','NDARTP382XC8','NDARLJ581GD7','NDARGS283RM9','NDARRED356WS', 'NDARHUG535MO'); %

num_conditions = 2;

%% Analysis Parameters
method = 'weight'; % 'choose' or 'weight'
mode = 'SPEECH VS NOISE'; % 'SPEECH', 'NOISE', or 'BOTH' (add 'NO BREATH' for no breath)
analysis_type = 'collapsed attend and masker PFC time constant';
area = 'PFC'; % 'PFC' or 'STG'
statistic_to_plot = 'mean'; % 'mean' or 'beta'


%% Block Averages

if user == 'Ben'
    GroupResults = readtable(append('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\Group Results SRM-NIRS-EEG-1 ',analysis_type,' ',mode, '.csv'),'Format','auto');
    addpath('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\errorbar_files\errorbar_files');
else
    GroupResults = readtable(append('/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/Group Results SRM-NIRS-EEG-1 ',analysis_type,' ',mode, '.csv'),'Format','auto');
    addpath('/home/ben/Documents/GitHub/SRM-NIRS-EEG/errorbar_files/errorbar_files');
end

conditions = unique(GroupResults.Condition);
conditions(string(conditions) == 'Exhale') = [];
conditions(string(conditions) == 'Inhale') = [];
conditions(string(conditions) == 'Hold') = [];

% channels
dlpfc_channels = [1,2,3,4,5,6];
stg_channels = 7:14;
left_stg_channels = 11:14;
right_stg_channels = 7:10;

% time
epoch_time_limits = [-5,20];

all_block_averages = [];
for isubject = 1:size(subject_ID,1)
    
    if user == 'Ben'
        this_subject_table = readtable('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\' + string(subject_ID(isubject,:)) + ' block averages ' + mode +'.csv');
    else
        this_subject_table = readtable('/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/' + string(subject_ID(isubject,:)) + ' block averages ' + mode +'.csv');
    end
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

%% Plot block averages
figure;
%ITD Conditions
ymin = -0.3;
ymax = 0.6;
ymin = -0.125;
ymax = 0.15;
time = linspace(epoch_time_limits(1),epoch_time_limits(2),size(all_block_averages,4));
plotting_fs = 5;
lineprop_list = {'-k',{'or','markerfacecolor',[1,1,1],'MarkerIndices',1:plotting_fs:length(time)},{'or','markerfacecolor',[1,0,0],'MarkerIndices',1:plotting_fs:length(time)},{'or','markerfacecolor',[1,1,1],'MarkerIndices',1:plotting_fs:length(time)},{'or','markerfacecolor',[1,0,0],'MarkerIndices',1:plotting_fs:length(time)}};
hold on
for icondition = 2:3
    this_lineprop = lineprop_list(icondition);
    shadedErrorBar(time,squeeze(nanmean(block_averages_to_plot(:,icondition,:),1)),squeeze(nanstd(block_averages_to_plot(:,icondition,:),[],1))./(sqrt(size(subject_ID,1))-1),'lineProps',this_lineprop{1,1});
    hold on
end

legend({'Noise','Speech'},'FontSize',14,'AutoUpdate','off')
ylim([ymin,ymax])
xline(0,'LineWidth',2)
xline(12.8,'LineWidth',2)

