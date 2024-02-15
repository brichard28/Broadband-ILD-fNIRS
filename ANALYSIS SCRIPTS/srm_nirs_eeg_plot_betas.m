%% srm_nirs_eeg_plot_betas.m
user = 'Laptop';
method = 'weight'; % 'choose' or 'weight'
analysis_type = 'collapsed attend and masker PFC time constant';
if user == 'Laptop'
    GroupResults = readtable(append('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\Group Results SRM-NIRS-EEG-1 ',analysis_type,' SPEECH.csv'),'Format','auto');
    addpath('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\errorbar_files\errorbar_files');
else
    GroupResults = readtable(append('/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/Group Results SRM-NIRS-EEG-1 ',analysis_type ,' SPEECH.csv'),'Format','auto');
    addpath('/home/ben/Documents/GitHub/SRM-NIRS-EEG/errorbar_files/errorbar_files');
end
subjects = unique(GroupResults.ID);
subjects(string(subjects) == 'NDARLM531OY3') = [];
subjects(string(subjects) == 'NDARGT639XS6') = [];

channels = unique(GroupResults.ch_name);
channels = string(channels);
channels(contains(channels,'hbr')) = [];

conditions = unique(GroupResults.Condition);
conditions(string(conditions) == 'Exhale') = [];
conditions(string(conditions) == 'Inhale') = [];
conditions(string(conditions) == 'Hold') = [];
if contains(analysis_type,'collapsed attend and masker')
    conditions([4 5]) = conditions([5 4]);
end
% channels
dlpfc_channels = [1,2,3,4,5,6];
stg_channels = 7:14;
left_stg_channels = 11:14;
right_stg_channels = 7:10;

% conditions
attend_left_conditions = [1,5,7,9,11,15,17,19]; % not control
attend_right_conditions = [3,6,8,10,13,16,18,20]; % not control
speech_conditions = [11,13,15:20]; % excluding control
noise_conditions = [1,3,5:10]; % excluding control
noise_control_conditions = [2,4]; % target left, target right
speech_control_conditions = [12,14]; % target left, target right

% time
epoch_time_limits = [-5,35];

all_betas = [];
all_ses = [];
for isubject = 1:length(subjects)
    for ichannel = 1:length(channels)
        for icondition = 1:length(conditions)
            all_betas(isubject,ichannel,icondition) = GroupResults.theta(string(GroupResults.ID) == string(subjects(isubject)) & string(GroupResults.ch_name) == string(channels(ichannel)) & string(GroupResults.Condition) == string(conditions(icondition)) & string(GroupResults.Chroma) == "hbo");
            all_ses(isubject,ichannel,icondition) = GroupResults.se(string(GroupResults.ID) == string(subjects(isubject)) & string(GroupResults.ch_name) == string(channels(ichannel)) & string(GroupResults.Condition) == string(conditions(icondition)) & string(GroupResults.Chroma) == "hbo");

        end
    end
end

all_stg_betas = all_betas(:,stg_channels,:);
all_pfc_betas = all_betas(:,dlpfc_channels,:);

all_stg_ses = all_ses(:,stg_channels,:);
all_pfc_ses = all_ses(:,dlpfc_channels,:);

all_block_averages = [];
for isubject = 1:length(subjects)
    if user == 'Laptop'
        this_subject_table = readtable('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\' + string(subjects(isubject)) + ' block averages.csv');
        this_epochs_deleted = readtable('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\' + string(subjects(isubject)) + ' epochs deleted.csv');
    else
        this_subject_table = readtable('/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/' + string(subjects(isubject)) + ' block averages.csv');
        this_epochs_deleted = readtable('/home/ben/Documents/GitHub/SRM-NIRS-EEG/RESULTS DATA/' + string(subjects(isubject)) + ' epochs deleted.csv');
    end
    if numel(this_epochs_deleted) > 0
        this_epochs_deleted = table2array(this_epochs_deleted(2:end,2));
    else
        this_epochs_deleted = [];
    end
    all_epochs = this_subject_table.epoch;
    for ichannel = 1:length(channels)
        for icondition = 1:length(conditions)
            this_condition_epochs = this_subject_table.epoch;
            this_condition_epochs = this_condition_epochs(string(this_subject_table.condition) == conditions(icondition));
            unique_epochs = unique(this_condition_epochs);
            curr_hbo = eval('this_subject_table.' + erase(channels(ichannel),' hbo') + 'Hbo');
            z = [];
            for i = 1:length(unique_epochs)
                if ~ismember(i,this_epochs_deleted)% is not rejected
                    z(i,:) = curr_hbo(all_epochs == unique_epochs(i));
                else
                    z(i,:) = nan*ones(121,1);
                end
            end
            all_block_averages(isubject,ichannel,icondition,:) = nanmean(z,1);
        end
    end
end

all_stg_block_averages = all_block_averages(:,stg_channels,:,:);
all_pfc_block_averages = all_block_averages(:,dlpfc_channels,:,:);

% Tak channel with maximum beta in control conditions
control_condition = find(string(conditions) == 'control');
[~,which_channel_stg] = max(squeeze(all_stg_betas(:,:,control_condition)),[],2);
[~,which_channel_pfc] = max(squeeze(all_pfc_betas(:,:,control_condition)),[],2);
if method == 'choose'
    % choose the single channel with the strongest response in the control
    % condition, for each area
    for isubject = 1:length(subjects)
        % Store beta values
        stg_betas_to_plot(isubject,:) = all_stg_betas(isubject,which_channel_stg(isubject),:);
        pfc_betas_to_plot(isubject,:) = all_pfc_betas(isubject,which_channel_pfc(isubject),:);
    
        % Store block averages
        stg_block_averages_to_plot(isubject,:,:) = squeeze(all_stg_block_averages(isubject,which_channel_stg(isubject),:,:));
        pfc_block_averages_to_plot(isubject,:,:) = squeeze(all_pfc_block_averages(isubject,which_channel_pfc(isubject),:,:));
        %stg_block_averages_to_plot(isubject,:) = ;
        %
    end
elseif method == 'weight'
    % weight each beta value by the inverse of the standard error of the
    % GLM fit, and include all channels in the analysis
    for isubject = 1:length(subjects)
        % Calculate beta values
        %stg_betas_to_plot(isubject,:) = mean(all_stg_betas(isubject,:,:),[1,2]);
        these_stg_betas = [];
        for ichannelstg = 1:size(all_stg_betas,2)
            these_stg_betas = cat(2,these_stg_betas, squeeze(all_stg_betas(isubject,ichannelstg,:).*(1/all_stg_ses(isubject,ichannelstg,:))));
        end
        stg_betas_to_plot(isubject,:) = mean(these_stg_betas,2);
        %pfc_betas_to_plot(isubject,:,:) = mean(all_pfc_betas(isubject,:,:),[1,2]);
        these_pfc_betas = [];
        for ichannelpfc = 1:size(all_pfc_betas,2)
            these_pfc_betas = cat(2,these_pfc_betas, squeeze(all_pfc_betas(isubject,ichannelpfc,:).*(1/all_pfc_ses(isubject,ichannelpfc,:))));
        end
        pfc_betas_to_plot(isubject,:) = mean(these_pfc_betas,2);

        % Store block averages
        stg_block_averages_to_plot(isubject,:,:) = squeeze(nanmean(all_stg_block_averages(isubject,:,:,:),2));
        pfc_block_averages_to_plot(isubject,:,:) = squeeze(nanmean(all_pfc_block_averages(isubject,:,:,:),2));

    end
end

% normalize
%pfc_betas_to_plot = (pfc_betas_to_plot - min(pfc_betas_to_plot,[],'all'))/(max(pfc_betas_to_plot,[],'all') - min(pfc_betas_to_plot,[],'all'));
%stg_betas_to_plot = (stg_betas_to_plot - min(stg_betas_to_plot,[],'all'))/(max(stg_betas_to_plot,[],'all') - min(stg_betas_to_plot,[],'all'));

if contains(analysis_type,'collapsed attend and masker') % compare across condition
    % STG
    figure;
    bar(squeeze(mean(stg_betas_to_plot(:,2:end),1)))
    hold on
    errorbar(1:2,squeeze(mean(stg_betas_to_plot(:,2:3),1)),squeeze(std(stg_betas_to_plot(:,2:3),[],1))./(sqrt(length(subjects))-1),'k','LineWidth',2)
    errorbar(3:4,squeeze(mean(stg_betas_to_plot(:,4:end),1)),squeeze(std(stg_betas_to_plot(:,4:end),[],1))./(sqrt(length(subjects))-1),'k','LineWidth',2)
    %plot(stg_betas_to_plot(:,2:end)','Color',[0.7 0.7 0.7])
    xlabel('Condition','FontSize',18)
    ylabel('Beta (AU)','FontSize',18)
    xticklabels({'Small ITD','Large ITD','Natural ILD','Broadband ILD'})
    title('STG','FontSize',18)
    ax = gca;
    ax.LineWidth = 2;
    ylim([-0.4, 1])
    % Statistics
    [p,tbl,stats] = anova1(stg_betas_to_plot(:,2:end));
    c = multcompare(stats,'Display','off');
    tbl = array2table(c,"VariableNames", ...
        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    p_values = c(:,6);
    significant_comparisons = find(p_values < 0.05);
    for i = 1:length(significant_comparisons)
        group_a = c(significant_comparisons(i),1);
        group_b = c(significant_comparisons(i),2);
        r = 0.05 + (0.1-0.05) .* rand(1,1);
        line([group_a,group_b],[max(squeeze(mean(stg_betas_to_plot(:,2:end),1)))+r,max(squeeze(mean(stg_betas_to_plot(:,2:end),1)))+r],'Color','k')
        text(mean([group_a,group_b]),max(squeeze(mean(stg_betas_to_plot(:,2:end),1)))+0.15,'*','FontSize',24)
    end

    % PFC
    figure;
    bar(squeeze(mean(pfc_betas_to_plot(:,2:end),1)))
    hold on
    errorbar(1:2,squeeze(mean(pfc_betas_to_plot(:,2:3),1)),squeeze(std(pfc_betas_to_plot(:,2:3),[],1))./(sqrt(length(subjects))-1),'k','LineWidth',2)
    errorbar(3:4,squeeze(mean(pfc_betas_to_plot(:,4:end),1)),squeeze(std(pfc_betas_to_plot(:,4:end),[],1))./(sqrt(length(subjects))-1),'k','LineWidth',2)

    %plot(pfc_betas_to_plot(:,2:end)','Color',[0.7 0.7 0.7])
    xlabel('Condition','FontSize',18)
    ylabel('Beta (AU)','FontSize',18)
    xticklabels({'Small ITD','Large ITD','Natural ILD','Broadband ILD'})
    title('PFC','FontSize',18)
    ax = gca;
    ax.LineWidth = 2;
    ylim([-0.4, 1])
    % Statistics
    [p,tbl,stats] = anova1(pfc_betas_to_plot(:,2:end));
    c = multcompare(stats,'Display','off');
    tbl = array2table(c,"VariableNames", ...
        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    p_values = c(:,6);
    significant_comparisons = find(p_values < 0.05);
    for i = 1:length(significant_comparisons)
        group_a = c(significant_comparisons(i),1);
        group_b = c(significant_comparisons(i),2);
        r = 0.05 + (0.1-0.05) .* rand(1,1);
        line([group_a,group_b],[max(squeeze(mean(pfc_betas_to_plot(:,2:end),1)))+r,max(squeeze(mean(pfc_betas_to_plot(:,2:end),1)))+r],'Color','k')
        text(mean([group_a,group_b]),max(squeeze(mean(pfc_betas_to_plot(:,2:end),1)))+0.15,'*','FontSize',24)
    end
elseif contains(analysis_type,'collapsed attend and condition')
    % STG
    figure;
    bar(squeeze(mean(stg_betas_to_plot(:,2:end),1)))
    hold on
    errorbar(1:2,squeeze(mean(stg_betas_to_plot(:,2:end),1)),squeeze(std(stg_betas_to_plot(:,2:end),[],1))./(sqrt(length(subjects))-1),'k','LineWidth',2)
    %plot(stg_betas_to_plot(:,2:end)','Color',[0.7 0.7 0.7])
    xlabel('Condition','FontSize',18)
    ylabel('Beta (AU)','FontSize',18)
    xticklabels({'Noise Masker','Speech Masker'})
    title('STG','FontSize',18)
    ax = gca;
    ax.LineWidth = 2;
    % Statistics
    [p,tbl,stats] = anova2(stg_betas_to_plot(:,2:end),1,'off');
    c = multcompare(stats,'Display','off');
    tbl = array2table(c,"VariableNames", ...
        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    p_values = c(:,6);
    significant_comparisons = find(p_values < 0.05);
    for i = 1:length(significant_comparisons)
        group_a = c(significant_comparisons(i),1);
        group_b = c(significant_comparisons(i),2);
        r = 0.05 + (0.1-0.05) .* rand(1,1);
        line([group_a,group_b],[max(squeeze(mean(stg_betas_to_plot(:,2:end),1)))+r,max(squeeze(mean(stg_betas_to_plot(:,2:end),1)))+r],'Color','k')
        text(mean([group_a,group_b]),max(squeeze(mean(stg_betas_to_plot(:,2:end),1)))+0.15,'*','FontSize',24)
    end

    % PFC
    figure;
    bar(squeeze(mean(pfc_betas_to_plot(:,2:end),1)))
    hold on
    errorbar(1:2,squeeze(mean(pfc_betas_to_plot(:,2:end),1)),squeeze(std(pfc_betas_to_plot(:,2:end),[],1))./(sqrt(length(subjects))-1),'k','LineWidth',2)
    %plot(pfc_betas_to_plot(:,2:end)','Color',[0.7 0.7 0.7])
    xlabel('Condition','FontSize',18)
    ylabel('Beta (AU)','FontSize',18)
    xticklabels({'Noise Masker','Speech Masker'})
    title('PFC','FontSize',18)
    ax = gca;
    ax.LineWidth = 2;
    % Statistics
    [p,tbl,stats] = anova2(pfc_betas_to_plot(:,2:end),1,'off');
    c = multcompare(stats,'Display','off');
    tbl = array2table(c,"VariableNames", ...
        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    p_values = c(:,6);
    significant_comparisons = find(p_values < 0.05);
    for i = 1:length(significant_comparisons)
        group_a = c(significant_comparisons(i),1);
        group_b = c(significant_comparisons(i),2);
        r = 0.05 + (0.1-0.05) .* rand(1,1);
        line([group_a,group_b],[max(squeeze(mean(pfc_betas_to_plot(:,2:end),1)))+r,max(squeeze(mean(pfc_betas_to_plot(:,2:end),1)))+r],'Color','k')
        text(mean([group_a,group_b]),max(squeeze(mean(pfc_betas_to_plot(:,2:end),1)))+0.15,'*','FontSize',24)
    end
end
%% Block Averages
if contains(analysis_type,'collapsed attend and masker') % compare across condition
    time = linspace(epoch_time_limits(1),epoch_time_limits(2),size(all_block_averages,4));
    lineprop_list = {'-k','-r','-g','-b','-m'};
    % STG
    figure;
    for icondition = 2:length(conditions)
        shadedErrorBar(time,squeeze(nanmean(stg_block_averages_to_plot(:,icondition,:),1)),squeeze(nanstd(stg_block_averages_to_plot(:,icondition,:),[],1))./(sqrt(length(subjects))-1),'lineProps',lineprop_list(icondition));
        hold on
    end
    xlabel('Time re stimulus onset (s)','FontSize',18)
    ylabel('\Delta HbO (AU)','FontSize',18)
    legend({'Small ITD','Large ITD','Natural ILD','Broadband ILD'},'AutoUpdate','off')
    title('STG Block Averages','FontSize',18)
    ylim([-4e5,10e5])
    xline(0,'LineWidth',1.5)
    xline(10,'LineWidth',1.5)
    % STG individual
    figure;
    for icondition = 2:length(conditions)
        for isubject = 1:length(subjects)
        plot(time,squeeze(stg_block_averages_to_plot(isubject,icondition,:)),string(lineprop_list(icondition)));
        hold on
        end
    end
    xlabel('Time re stimulus onset (s)','FontSize',18)
    ylabel('\Delta HbO (AU)','FontSize',18)
    %legend({'ITD50','ITD500','ILD10','ILD70n'})
    title('Individual STG Block Averages','FontSize',18)
    ylim([-4e5,10e5])

    % PFC
    figure;
    for icondition = 2:length(conditions)
%         if icondition == 2
%             shadedErrorBar(time,squeeze(nanmean(0.3*pfc_block_averages_to_plot(:,icondition,:),1)),squeeze(nanstd(0.3*pfc_block_averages_to_plot(:,icondition,:),[],1))./(sqrt(length(subjects))-1),'lineProps',lineprop_list(icondition));
% 
%         else
            shadedErrorBar(time,squeeze(nanmean(pfc_block_averages_to_plot(:,icondition,:),1)),squeeze(nanstd(pfc_block_averages_to_plot(:,icondition,:),[],1))./(sqrt(length(subjects))-1),'lineProps',lineprop_list(icondition));

%         end
        hold on
    end
    xlabel('Time re stimulus onset (s)','FontSize',18)
    ylabel('\Delta HbO (AU)','FontSize',18)
    legend({'Small ITD','Large ITD','Natural ILD','Broadband ILD'},'AutoUpdate','off')
    title('PFC Block Averages','FontSize',18)
    ylim([-4e5,10e5])
    xline(0,'LineWidth',1.5)
    xline(10,'LineWidth',1.5)
    % PFC individual
    figure;
    for icondition = 2:length(conditions)
        for isubject = 1:length(subjects)
        plot(time,squeeze(pfc_block_averages_to_plot(isubject,icondition,:)),string(lineprop_list(icondition)));
        hold on
        end
    end
    xlabel('Time re stimulus onset (s)','FontSize',18)
    ylabel('\Delta HbO (AU)','FontSize',18)
    %legend({'ITD50','ITD500','ILD10','ILD70n'})
    title('Individual PFC Block Averages','FontSize',18)
    ylim([-4e5,10e5])

end

%% Area Under the Curve
if contains(analysis_type,'collapsed attend and masker') % compare across condition
    for icondition = 2:length(conditions)
        for isubject = 1:length(subjects)
            % calculate AUC for all blocks
            [~,time_index_0] = min(abs(time - 0));
            [~,time_index_10] = min(abs(time - 10.8));
            this_response_STG = stg_block_averages_to_plot(isubject,icondition,time_index_0:time_index_10);
            this_response_PFC = pfc_block_averages_to_plot(isubject,icondition,time_index_0:time_index_10);
            x = time(time_index_0:time_index_10);
            all_AUC_stg(isubject,icondition) = trapz(x, this_response_STG);
            all_AUC_pfc(isubject,icondition) = trapz(x, this_response_PFC);
        end
    end
    figure;

    bar(squeeze(mean(all_AUC_stg(:,2:end),1)))
    hold on
    errorbar(1:2,squeeze(mean(all_AUC_stg(:,2:3),1)),squeeze(std(all_AUC_stg(:,2:3),[],1))./(sqrt(length(subjects))-1),'k','LineWidth',2)
    errorbar(3:4,squeeze(mean(all_AUC_stg(:,4:end),1)),squeeze(std(all_AUC_stg(:,4:end),[],1))./(sqrt(length(subjects))-1),'k','LineWidth',2)
    xlabel('Condition','FontSize',18)
    ylabel('AUC (AU)','FontSize',18)
    xticklabels({'Small ITD','Large ITD','Natural ILD','Broadband ILD'})
    title('Area Under the Curve STG')
    figure;
    bar(squeeze(mean(all_AUC_pfc(:,2:end),1)))
    hold on
    errorbar(1:2,squeeze(mean(all_AUC_pfc(:,2:3),1)),squeeze(std(all_AUC_pfc(:,2:3),[],1))./(sqrt(length(subjects))-1),'k','LineWidth',2)
    errorbar(3:4,squeeze(mean(all_AUC_pfc(:,4:end),1)),squeeze(std(all_AUC_pfc(:,4:end),[],1))./(sqrt(length(subjects))-1),'k','LineWidth',2)
    xlabel('Condition','FontSize',18)
    ylabel('AUC (AU)','FontSize',18)
    xticklabels({'Small ITD','Large ITD','Natural ILD','Broadband ILD'})
    title('Area Under the Curve PFC')

end


% save
save(['C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\',analysis_type,' SPEECH results.mat'],'stg_betas_to_plot','pfc_betas_to_plot','all_AUC_stg','all_AUC_pfc');


%% STG vs. PFC


%% OLD CODE
% % plot speech attend left vs. speech attend right
% all_stg_noise_betas = all_betas(:,stg_channels,4:9);
% all_stg_speech_betas = all_betas(:,stg_channels,10:15);


%
% %% INCLUDING ALL CHANNELS
% % Speech plot Left Hemisphere (attend left vs. attend right)
% upper_ylim = 1;
% lower_ylim = -1;
% figure;
% subplot(2,2,1)
% attend_left_betas = all_betas(:,left_stg_channels,intersect(speech_conditions,attend_left_conditions));
% attend_right_betas = all_betas(:,left_stg_channels,intersect(speech_conditions,attend_right_conditions));
% x_axis_attend_left = 1:3:10;
% x_axis_attend_right = 2:3:11;
% p1 = violinplot_T(reshape(attend_left_betas,[size(attend_left_betas,1)*size(attend_left_betas,2),size(attend_left_betas,3)]),{'1','2','3','4'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% hold on;p2 =  violinplot_T(reshape(attend_right_betas,[size(attend_right_betas,1)*size(attend_right_betas,2),size(attend_right_betas,3)]),{'1','2','3','4'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
% legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
% xticks(1.5:3:16.5)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Speech Masker Left STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
%
%
% % Speech plot Right Hemisphere (attend_left vs. attend_right)
% subplot(2,2,2)
% attend_left_betas = all_betas(:,right_stg_channels,intersect(speech_conditions,attend_left_conditions));
% attend_right_betas = all_betas(:,right_stg_channels,intersect(speech_conditions,attend_right_conditions));
% p1 = violinplot_T(reshape(attend_left_betas,[size(attend_left_betas,1)*size(attend_left_betas,2),size(attend_left_betas,3)]),{'1','2','3','4'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% hold on;p2 =  violinplot_T(reshape(attend_right_betas,[size(attend_right_betas,1)*size(attend_right_betas,2),size(attend_right_betas,3)]),{'1','2','3','4'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
% legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
% xticks(1.5:3:16.5)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Speech Masker Right STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
%
% % Noise plot Left Hemisphere (attend left vs. attend right)
% subplot(2,2,3)
% attend_left_betas = all_betas(:,left_stg_channels,intersect(noise_conditions,attend_left_conditions));
% attend_right_betas = all_betas(:,left_stg_channels,intersect(noise_conditions,attend_right_conditions));
% p1 = violinplot_T(reshape(attend_left_betas,[size(attend_left_betas,1)*size(attend_left_betas,2),size(attend_left_betas,3)]),{'1','2','3','4'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% hold on;p2 =  violinplot_T(reshape(attend_right_betas,[size(attend_right_betas,1)*size(attend_right_betas,2),size(attend_right_betas,3)]),{'1','2','3','4'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
% legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
% xticks(1.5:3:16.5)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Noise Masker Left STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
%
% % Noise plot Right hemisphere (attend left vs. attend right)
% subplot(2,2,4)
% attend_left_betas = all_betas(:,right_stg_channels,intersect(noise_conditions,attend_left_conditions));
% attend_right_betas = all_betas(:,right_stg_channels,intersect(noise_conditions,attend_right_conditions));
% p1 = violinplot_T(reshape(attend_left_betas,[size(attend_left_betas,1)*size(attend_left_betas,2),size(attend_left_betas,3)]),{'1','2','3','4'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% hold on;p2 =  violinplot_T(reshape(attend_right_betas,[size(attend_right_betas,1)*size(attend_right_betas,2),size(attend_right_betas,3)]),{'1','2','3','4'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
% legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
% xticks(1.5:3:16.5)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Noise Masker Right STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
%
% sgtitle('REGULAR DONT TOUCH ALL CHANNELS REPRESENTED HERE')
%
%
%
%
%
%
% %% CHOOSING CHANNEL ON EACH SIDE BASED ON ATTEND LEFT or ATTEND RIGHT ITD 500 SPEECH NOISE CONTRAST
%
% left_stg_control_conditions = all_betas(:,left_stg_channels,speech_control_conditions(1)) - all_betas(:,left_stg_channels,noise_control_conditions(1));
% right_stg_control_conditions = all_betas(:,right_stg_channels,speech_control_conditions(1)) - all_betas(:,right_stg_channels,noise_control_conditions(1));
%
% [~,channels_to_choose_by_subject_left_stg] = max(left_stg_control_conditions,[],2);
% [~,channels_to_choose_by_subject_right_stg] = max(right_stg_control_conditions,[],2);
%
% % Speech plot Left Hemisphere (attend left vs. attend right)
% figure;
% subplot(2,2,1)
% attend_left_betas = [];
% attend_right_betas = [];
% for isubject = 1:length(subjects)
%     attend_left_betas(isubject,:,:) = all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg(isubject)),intersect(speech_conditions,attend_left_conditions));
%     attend_right_betas(isubject,:,:) = all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg(isubject)),intersect(speech_conditions,attend_right_conditions));
% end
% x_axis_attend_left = 1:3:10;
% x_axis_attend_right = 2:3:11;
% p1 = violinplot_T(attend_left_betas,{'1','2','3'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% hold on;p2 =  violinplot_T(attend_right_betas,{'1','2','3'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
% legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
% xticks(1.5:3:16.5)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Speech Masker Left STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
%
% % Speech plot Right Hemisphere (attend_left vs. attend_right)
% subplot(2,2,2)
% attend_left_betas = [];
% attend_right_betas = [];
% for isubject = 1:length(subjects)
%
% attend_left_betas(isubject,:,:) = all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg(isubject)),intersect(speech_conditions,attend_left_conditions));
% attend_right_betas(isubject,:,:) = all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg(isubject)),intersect(speech_conditions,attend_right_conditions));
% end
% p1 = violinplot_T(attend_left_betas,{'1','2','3','4'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% hold on;p2 =  violinplot_T(attend_right_betas,{'1','2','3','4'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
% legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
% xticks(1.5:3:16.5)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Speech Masker Right STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
% % Noise plot Left Hemisphere (attend left vs. attend right)
% subplot(2,2,3)
% attend_left_betas = [];
% attend_right_betas = [];
% for isubject = 1:length(subjects)
%
% attend_left_betas(isubject,:,:) = all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg(isubject)),intersect(noise_conditions,attend_left_conditions));
% attend_right_betas(isubject,:,:) = all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg(isubject)),intersect(noise_conditions,attend_right_conditions));
% end
% p1 = violinplot_T(attend_left_betas,{'1','2','3','4'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% hold on;p2 =  violinplot_T(attend_right_betas,{'1','2','3','4'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
% legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
% xticks(1.5:3:16.5)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Noise Masker Left STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
% % Noise plot Right hemisphere (attend left vs. attend right)
% subplot(2,2,4)
% attend_left_betas = [];
% attend_right_betas = [];
% for isubject = 1:length(subjects)
%
% attend_left_betas(isubject,:,:) = all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg(isubject)),intersect(noise_conditions,attend_left_conditions));
% attend_right_betas(isubject,:,:) = all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg(isubject)),intersect(noise_conditions,attend_right_conditions));
% end
% p1 = violinplot_T(attend_left_betas,{'1','2','3','4'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% hold on;p2 =  violinplot_T(attend_right_betas,{'1','2','3','4'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
% legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
% xticks(1.5:3:16.5)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Noise Masker Right STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
%
% sgtitle('Channel Chosen from ITD500 Control Condition Speech vs. Noise ATTEND LEFT')
%
% %% COLLAPSED ACROSS LEFT AND RIGHT ATTEND LEFT CHANNEL CHOICE STG
% % Left Hemisphere Speech
% figure;
% subplot(2,2,1)
% betas_to_plot = [];
% for isubject = 1:length(subjects)
%     a =  all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg(isubject)),intersect(speech_conditions,attend_left_conditions));
%     b = all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg(isubject)),intersect(speech_conditions,attend_right_conditions));
%     betas_to_plot(isubject,:,:) = mean([a,b],2);
% end
% p1 = violinplot_T(betas_to_plot,{'1','2','3','4'},[1,2,3,4],'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% xticks(1:4)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Speech Masker Left STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
% % Right Hemisphere Speech
% subplot(2,2,2)
% betas_to_plot = [];
% for isubject = 1:length(subjects)
%     a =  all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg(isubject)),intersect(speech_conditions,attend_left_conditions));
%     b = all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg(isubject)),intersect(speech_conditions,attend_right_conditions));
%     betas_to_plot(isubject,:,:) = mean([a,b],2);
% end
% p1 = violinplot_T(betas_to_plot,{'1','2','3','4'},[1,2,3,4],'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% xticks(1:4)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Speech Masker Left STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
% % Left Hemisphere Noise
% subplot(2,2,3)
% betas_to_plot = [];
% for isubject = 1:length(subjects)
%     a =  all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg(isubject)),intersect(noise_conditions,attend_left_conditions));
%     b = all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg(isubject)),intersect(noise_conditions,attend_right_conditions));
%     betas_to_plot(isubject,:,:) = mean([a,b],2);
% end
% p1 = violinplot_T(betas_to_plot,{'1','2','3','4'},[1,2,3,4],'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% xticks(1:4)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Noise Masker Left STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
% % Right Hemisphere Noise
% subplot(2,2,4)
% betas_to_plot = [];
% for isubject = 1:length(subjects)
%     a =  all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg(isubject)),intersect(noise_conditions,attend_left_conditions));
%     b = all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg(isubject)),intersect(noise_conditions,attend_right_conditions));
%     betas_to_plot(isubject,:,:) = mean([a,b],2);
% end
% p1 = violinplot_T(betas_to_plot,{'1','2','3','4'},[1,2,3,4],'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% xticks(1:4)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Noise Masker Right STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
% sgtitle('ATTEND RIGHT Channel Choice, collapsed attend left and right')
%
% %% CHOOSING CHANNEL ON EACH SIDE BASED ON ATTEND LEFT or ATTEND RIGHT ITD 500 SPEECH NOISE CONTRAST
%
% left_stg_control_conditions = all_betas(:,left_stg_channels,speech_control_conditions(2)) - all_betas(:,left_stg_channels,noise_control_conditions(2));
% right_stg_control_conditions = all_betas(:,right_stg_channels,speech_control_conditions(2)) - all_betas(:,right_stg_channels,noise_control_conditions(2));
%
% [~,channels_to_choose_by_subject_left_stg] = max(left_stg_control_conditions,[],2);
% [~,channels_to_choose_by_subject_right_stg] = max(right_stg_control_conditions,[],2);
%
% % Speech plot Left Hemisphere (attend left vs. attend right)
% figure;
% subplot(2,2,1)
% attend_left_betas = [];
% attend_right_betas = [];
% for isubject = 1:length(subjects)
%     attend_left_betas(isubject,:,:) = all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg(isubject)),intersect(speech_conditions,attend_left_conditions));
%     attend_right_betas(isubject,:,:) = all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg(isubject)),intersect(speech_conditions,attend_right_conditions));
% end
% x_axis_attend_left = 1:3:10;
% x_axis_attend_right = 2:3:11;
% p1 = violinplot_T(attend_left_betas,{'1','2','3'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% hold on;p2 =  violinplot_T(attend_right_betas,{'1','2','3'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
% legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
% xticks(1.5:3:16.5)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Speech Masker Left STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
%
% % Speech plot Right Hemisphere (attend_left vs. attend_right)
% subplot(2,2,2)
% attend_left_betas = [];
% attend_right_betas = [];
% for isubject = 1:length(subjects)
%
% attend_left_betas(isubject,:,:) = all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg(isubject)),intersect(speech_conditions,attend_left_conditions));
% attend_right_betas(isubject,:,:) = all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg(isubject)),intersect(speech_conditions,attend_right_conditions));
% end
% p1 = violinplot_T(attend_left_betas,{'1','2','3','4'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% hold on;p2 =  violinplot_T(attend_right_betas,{'1','2','3','4'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
% legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
% xticks(1.5:3:16.5)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Speech Masker Right STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
% % Noise plot Left Hemisphere (attend left vs. attend right)
% subplot(2,2,3)
% attend_left_betas = [];
% attend_right_betas = [];
% for isubject = 1:length(subjects)
%
% attend_left_betas(isubject,:,:) = all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg(isubject)),intersect(noise_conditions,attend_left_conditions));
% attend_right_betas(isubject,:,:) = all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg(isubject)),intersect(noise_conditions,attend_right_conditions));
% end
% p1 = violinplot_T(attend_left_betas,{'1','2','3','4'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% hold on;p2 =  violinplot_T(attend_right_betas,{'1','2','3','4'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
% legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
% xticks(1.5:3:16.5)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Noise Masker Left STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
% % Noise plot Right hemisphere (attend left vs. attend right)
% subplot(2,2,4)
% attend_left_betas = [];
% attend_right_betas = [];
% for isubject = 1:length(subjects)
%
% attend_left_betas(isubject,:,:) = all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg(isubject)),intersect(noise_conditions,attend_left_conditions));
% attend_right_betas(isubject,:,:) = all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg(isubject)),intersect(noise_conditions,attend_right_conditions));
% end
% p1 = violinplot_T(attend_left_betas,{'1','2','3','4'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% hold on;p2 =  violinplot_T(attend_right_betas,{'1','2','3','4'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
% legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
% xticks(1.5:3:16.5)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Noise Masker Right STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
%
% sgtitle('Channel Chosen from ITD500 Control Condition Speech vs. Noise ATTEND RIGHT')
%
% %% COLLAPSED ACROSS LEFT AND RIGHT ATTEND RIGHT CHANNEL CHOICE STG
% % Left Hemisphere Speech
% figure;
% subplot(2,2,1)
% betas_to_plot = [];
% for isubject = 1:length(subjects)
%     a =  all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg(isubject)),intersect(speech_conditions,attend_left_conditions));
%     b = all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg(isubject)),intersect(speech_conditions,attend_right_conditions));
%     betas_to_plot(isubject,:,:) = mean([a,b],2);
% end
% p1 = violinplot_T(betas_to_plot,{'1','2','3','4'},[1,2,3,4],'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% xticks(1:4)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Speech Masker Left STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
% % Right Hemisphere Speech
% subplot(2,2,2)
% betas_to_plot = [];
% for isubject = 1:length(subjects)
%     a =  all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg(isubject)),intersect(speech_conditions,attend_left_conditions));
%     b = all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg(isubject)),intersect(speech_conditions,attend_right_conditions));
%     betas_to_plot(isubject,:,:) = mean([a,b],2);
% end
% p1 = violinplot_T(betas_to_plot,{'1','2','3','4'},[1,2,3,4],'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% xticks(1:4)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Speech Masker Left STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
% % Left Hemisphere Noise
% subplot(2,2,3)
% betas_to_plot = [];
% for isubject = 1:length(subjects)
%     a =  all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg(isubject)),intersect(noise_conditions,attend_left_conditions));
%     b = all_betas(isubject,left_stg_channels(channels_to_choose_by_subject_left_stg(isubject)),intersect(noise_conditions,attend_right_conditions));
%     betas_to_plot(isubject,:,:) = mean([a,b],2);
% end
% p1 = violinplot_T(betas_to_plot,{'1','2','3','4'},[1,2,3,4],'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% xticks(1:4)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Noise Masker Left STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
% % Right Hemisphere Noise
% subplot(2,2,4)
% betas_to_plot = [];
% for isubject = 1:length(subjects)
%     a =  all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg(isubject)),intersect(noise_conditions,attend_left_conditions));
%     b = all_betas(isubject,right_stg_channels(channels_to_choose_by_subject_right_stg(isubject)),intersect(noise_conditions,attend_right_conditions));
%     betas_to_plot(isubject,:,:) = mean([a,b],2);
% end
% p1 = violinplot_T(betas_to_plot,{'1','2','3','4'},[1,2,3,4],'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% xticks(1:4)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Noise Masker Right STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
% sgtitle('ATTEND RIGHT Channel Choice, collapsed attend left and right')
%
%
% %% Prefrontal Cortex
% %% INCLUDING ALL CHANNELS
% % Speech plot (attend left vs. attend right)
% upper_ylim = 0.3;
% lower_ylim = -0.3;
% figure;
% subplot(2,1,1)
% attend_left_betas = all_betas(:,dlpfc_channels,intersect(speech_conditions,attend_left_conditions));
% attend_right_betas = all_betas(:,dlpfc_channels,intersect(speech_conditions,attend_right_conditions));
% x_axis_attend_left = 1:3:10;
% x_axis_attend_right = 2:3:11;
% p1 = violinplot_T(reshape(attend_left_betas,[size(attend_left_betas,1)*size(attend_left_betas,2),size(attend_left_betas,3)]),{'1','2','3','4'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% hold on;p2 =  violinplot_T(reshape(attend_right_betas,[size(attend_right_betas,1)*size(attend_right_betas,2),size(attend_right_betas,3)]),{'1','2','3','4'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
% legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
% xticks(1.5:3:16.5)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Speech Masker','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
%
%
% % Noise plot (attend_left vs. attend_right)
% subplot(2,1,2)
% attend_left_betas = all_betas(:,dlpfc_channels,intersect(noise_conditions,attend_left_conditions));
% attend_right_betas = all_betas(:,dlpfc_channels,intersect(noise_conditions,attend_right_conditions));
% p1 = violinplot_T(reshape(attend_left_betas,[size(attend_left_betas,1)*size(attend_left_betas,2),size(attend_left_betas,3)]),{'1','2','3','4'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% hold on;p2 =  violinplot_T(reshape(attend_right_betas,[size(attend_right_betas,1)*size(attend_right_betas,2),size(attend_right_betas,3)]),{'1','2','3','4'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
% legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
% xticks(1.5:3:16.5)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Noise Masker','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
%
% sgtitle('PFC All Channels')
%
%
% %% CHOOSE CHANNELS BASED ON ITD500 CONTROL CONDITIONS
% dlpfc_control_conditions = all_betas(:,dlpfc_channels,speech_control_conditions(2)) - all_betas(:,dlpfc_channels,noise_control_conditions(2));
%
% [~,channels_to_choose_by_subject_pfc] = max(dlpfc_control_conditions,[],2);
%
% % Speech plot (attend left vs. attend right)
% upper_ylim = 0.3;
% lower_ylim = -0.3;
% figure;
% subplot(2,1,1)
% attend_left_betas = [];
% attend_right_betas = [];
% for isubject = 1:length(subjects)
% attend_left_betas(isubject,:,:) = all_betas(isubject,dlpfc_channels(channels_to_choose_by_subject_pfc(isubject)),intersect(speech_conditions,attend_left_conditions));
% attend_right_betas(isubject,:,:) = all_betas(isubject,dlpfc_channels(channels_to_choose_by_subject_pfc(isubject)),intersect(speech_conditions,attend_right_conditions));
% end
% x_axis_attend_left = 1:3:10;
% x_axis_attend_right = 2:3:11;
% p1 = violinplot_T(attend_left_betas,{'1','2','3','4'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% hold on;p2 =  violinplot_T(attend_right_betas,{'1','2','3','4'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
% legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
% xticks(1.5:3:16.5)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Speech Masker','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
%
%
% % Noise plot (attend_left vs. attend_right)
% subplot(2,1,2)
% attend_left_betas = [];
% attend_right_betas = [];
% for isubject = 1:length(subjects)
% attend_left_betas(isubject,:,:) = all_betas(isubject,dlpfc_channels(channels_to_choose_by_subject_pfc(isubject)),intersect(noise_conditions,attend_left_conditions));
% attend_right_betas(isubject,:,:) = all_betas(isubject,dlpfc_channels(channels_to_choose_by_subject_pfc(isubject)),intersect(noise_conditions,attend_right_conditions));
% end
% p1 = violinplot_T(attend_left_betas,{'1','2','3','4'},x_axis_attend_left,'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% hold on;p2 =  violinplot_T(attend_right_betas,{'1','2','3','4'},x_axis_attend_right,'ViolinColor',[1 0 0],'ViolinAlpha',{0.5 0.02}');
% legend([p1(1,1).ViolinPlot,p2(1,1).ViolinPlot],{'Attend Left','Attend Right'})
% xticks(1.5:3:16.5)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Noise Masker','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
% sgtitle('PFC Channel Choice')
%
%
% %% COLLAPSED ACROSS LEFT AND RIGHT CHANNEL CHOICE PFC
% %% COLLAPSED ACROSS LEFT AND RIGHT ATTEND RIGHT CHANNEL CHOICE STG
% % Left Hemisphere Speech
% figure;
% subplot(2,2,1)
% betas_to_plot = [];
% for isubject = 1:length(subjects)
%     a =  all_betas(isubject,dlpfc_channels(channels_to_choose_by_subject_pfc(isubject)),intersect(speech_conditions,attend_left_conditions));
%     b = all_betas(isubject,dlpfc_channels(channels_to_choose_by_subject_pfc(isubject)),intersect(speech_conditions,attend_right_conditions));
%     betas_to_plot(isubject,:,:) = mean([a,b],2);
% end
% p1 = violinplot_T(betas_to_plot,{'1','2','3','4'},[1,2,3,4],'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% xticks(1:4)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Speech Masker Left STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
% % Right Hemisphere Speech
% subplot(2,2,2)
% betas_to_plot = [];
% for isubject = 1:length(subjects)
%     a =  all_betas(isubject,dlpfc_channels(channels_to_choose_by_subject_pfc(isubject)),intersect(speech_conditions,attend_left_conditions));
%     b = all_betas(isubject,dlpfc_channels(channels_to_choose_by_subject_pfc(isubject)),intersect(speech_conditions,attend_right_conditions));
%     betas_to_plot(isubject,:,:) = mean([a,b],2);
% end
% p1 = violinplot_T(betas_to_plot,{'1','2','3','4'},[1,2,3,4],'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% xticks(1:4)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Speech Masker Left STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
% % Left Hemisphere Noise
% subplot(2,2,3)
% betas_to_plot = [];
% for isubject = 1:length(subjects)
%     a =  all_betas(isubject,dlpfc_channels(channels_to_choose_by_subject_pfc(isubject)),intersect(noise_conditions,attend_left_conditions));
%     b = all_betas(isubject,dlpfc_channels(channels_to_choose_by_subject_pfc(isubject)),intersect(noise_conditions,attend_right_conditions));
%     betas_to_plot(isubject,:,:) = mean([a,b],2);
% end
% p1 = violinplot_T(betas_to_plot,{'1','2','3','4'},[1,2,3,4],'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% xticks(1:4)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Noise Masker Left STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
% % Right Hemisphere Noise
% subplot(2,2,4)
% betas_to_plot = [];
% for isubject = 1:length(subjects)
%     a =  all_betas(isubject,dlpfc_channels(channels_to_choose_by_subject_pfc(isubject)),intersect(noise_conditions,attend_left_conditions));
%     b = all_betas(isubject,dlpfc_channels(channels_to_choose_by_subject_pfc(isubject)),intersect(noise_conditions,attend_right_conditions));
%     betas_to_plot(isubject,:,:) = mean([a,b],2);
% end
% p1 = violinplot_T(betas_to_plot,{'1','2','3','4'},[1,2,3,4],'ViolinColor',[0 0 1],'ViolinAlpha',{0.5 0.02}');
% xticks(1:4)
% xticklabels({'ITD 500','ITD 50','ILD 10','ILD Nat'})
% title('Noise Masker Right STG','FontSize',18)
% xlabel('Condition','FontSize',18)
% ylabel('Beta','FontSize',18)
% ylim([lower_ylim upper_ylim])
%
% sgtitle('ATTEND RIGHT Channel Choice, collapsed attend left and right')