%% srm_nirs_eeg_plot_block_averages

% Plot block averages for paper
speech_masker_data = readtable("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_uncorr_block_average_speech_masker.csv",'format','auto');
noise_masker_data = readtable("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_uncorr_block_average_noise_masker.csv",'format','auto');
speech_masker_data = speech_masker_data(2:end,:);
noise_masker_data = noise_masker_data(2:end,:);

speech_masker_data.Properties.VariableNames  = {'S','Channel','TimeIndex','ITD50', 'ITD500', 'ILD70n', 'ILD10'};
noise_masker_data.Properties.VariableNames  = {'S','Channel','TimeIndex','ITD50', 'ITD500', 'ILD70n', 'ILD10'};

speech_masker_plot_data_pfc = [];
speech_masker_plot_data_stg = [];
noise_masker_plot_data_pfc = [];
noise_masker_plot_data_stg = [];

num_subjects = 30;

% PFC
for isubject = 0:(num_subjects -1)
    for ichannel = 0:5
        this_channel_speech_data = speech_masker_data(speech_masker_data.S == string(isubject) & speech_masker_data.Channel == string(ichannel),:);
        this_channel_noise_data = noise_masker_data(noise_masker_data.S == string(isubject) & noise_masker_data.Channel == string(ichannel),:);

        speech_masker_plot_data_pfc(isubject+1,ichannel+1,1,:) = this_channel_speech_data.ITD50;
        speech_masker_plot_data_pfc(isubject+1,ichannel+1,2,:) = this_channel_speech_data.ITD500;
        speech_masker_plot_data_pfc(isubject+1,ichannel+1,3,:) = this_channel_speech_data.ILD70n;
        speech_masker_plot_data_pfc(isubject+1,ichannel+1,4,:) = this_channel_speech_data.ILD10;

        noise_masker_plot_data_pfc(isubject+1,ichannel+1,1,:) = this_channel_noise_data.ITD50;
        noise_masker_plot_data_pfc(isubject+1,ichannel+1,2,:) = this_channel_noise_data.ITD500;
        noise_masker_plot_data_pfc(isubject+1,ichannel+1,3,:) = this_channel_noise_data.ILD70n;
        noise_masker_plot_data_pfc(isubject+1,ichannel+1,4,:) = this_channel_noise_data.ILD10;

    end

end

% STG

for isubject = 0:(num_subjects -1)
    for ichannel = 6:13
        this_channel_speech_data = speech_masker_data(speech_masker_data.S == string(isubject) & speech_masker_data.Channel == string(ichannel),:);
        this_channel_noise_data = noise_masker_data(noise_masker_data.S == string(isubject) & noise_masker_data.Channel == string(ichannel),:);

        speech_masker_plot_data_stg(isubject+1,ichannel-5,1,:) = this_channel_speech_data.ITD50;
        speech_masker_plot_data_stg(isubject+1,ichannel-5,2,:) = this_channel_speech_data.ITD500;
        speech_masker_plot_data_stg(isubject+1,ichannel-5,3,:) = this_channel_speech_data.ILD70n;
        speech_masker_plot_data_stg(isubject+1,ichannel-5,4,:) = this_channel_speech_data.ILD10;

        noise_masker_plot_data_stg(isubject+1,ichannel-5,1,:) = this_channel_noise_data.ITD50;
        noise_masker_plot_data_stg(isubject+1,ichannel-5,2,:) = this_channel_noise_data.ITD500;
        noise_masker_plot_data_stg(isubject+1,ichannel-5,3,:) = this_channel_noise_data.ILD70n;
        noise_masker_plot_data_stg(isubject+1,ichannel-5,4,:) = this_channel_noise_data.ILD10;

    end

end


% Figure prep stuff
epoch_start_time = -5;
epoch_end_time = 20;
num_timepoints = 255;
time = linspace(epoch_start_time,epoch_end_time,num_timepoints);
colors = {'r','b','m','g'};
ymin = -0.1;
ymax = 0.15;

figure;
tiledlayout(2,4)

% PFC, Speech masker
nexttile
hold on
shadedErrorBar(time,squeeze(mean(speech_masker_plot_data_pfc(:,:,1,:),[1,2,3],'omitnan')),squeeze(std(speech_masker_plot_data_pfc(:,:,1,:),[],[1,2,3],'omitnan'))./sqrt(num_subjects-1), 'lineProps',colors(1) )
shadedErrorBar(time,squeeze(mean(speech_masker_plot_data_pfc(:,:,2,:),[1,2,3],'omitnan')),squeeze(std(speech_masker_plot_data_pfc(:,:,2,:),[],[1,2,3],'omitnan'))./sqrt(num_subjects-1), 'lineProps',colors(2)  )
ylim([ymin,ymax])
xlim([epoch_start_time,epoch_end_time])
set(gca,'fontsize',14)


nexttile
hold on
shadedErrorBar(time,squeeze(mean(speech_masker_plot_data_pfc(:,:,3,:),[1,2,3],'omitnan')),squeeze(std(speech_masker_plot_data_pfc(:,:,3,:),[],[1,2,3],'omitnan'))./sqrt(num_subjects-1), 'lineProps',colors(3)  )
shadedErrorBar(time,squeeze(mean(speech_masker_plot_data_pfc(:,:,4,:),[1,2,3],'omitnan')),squeeze(std(speech_masker_plot_data_pfc(:,:,4,:),[],[1,2,3],'omitnan'))./sqrt(num_subjects-1), 'lineProps',colors(4)  )
ylim([ymin,ymax])
xlim([epoch_start_time,epoch_end_time])
set(gca,'fontsize',14)


% STG, Speech Masker
nexttile
hold on
shadedErrorBar(time,squeeze(mean(speech_masker_plot_data_stg(:,:,1,:),[1,2,3],'omitnan')),squeeze(std(speech_masker_plot_data_stg(:,:,1,:),[],[1,2,3],'omitnan'))./sqrt(num_subjects-1), 'lineProps',colors(1) )
shadedErrorBar(time,squeeze(mean(speech_masker_plot_data_stg(:,:,2,:),[1,2,3],'omitnan')),squeeze(std(speech_masker_plot_data_stg(:,:,2,:),[],[1,2,3],'omitnan'))./sqrt(num_subjects-1), 'lineProps',colors(2)  )
ylim([ymin,ymax])
xlim([epoch_start_time,epoch_end_time])
set(gca,'fontsize',14)


nexttile
hold on
shadedErrorBar(time,squeeze(mean(speech_masker_plot_data_stg(:,:,3,:),[1,2,3],'omitnan')),squeeze(std(speech_masker_plot_data_stg(:,:,3,:),[],[1,2,3],'omitnan'))./sqrt(num_subjects-1), 'lineProps',colors(3)  )
shadedErrorBar(time,squeeze(mean(speech_masker_plot_data_stg(:,:,4,:),[1,2,3],'omitnan')),squeeze(std(speech_masker_plot_data_stg(:,:,4,:),[],[1,2,3],'omitnan'))./sqrt(num_subjects-1), 'lineProps',colors(4)  )
ylim([ymin,ymax])
xlim([epoch_start_time,epoch_end_time])
set(gca,'fontsize',14)


% PFC, noise masker
nexttile
hold on
shadedErrorBar(time,squeeze(mean(noise_masker_plot_data_pfc(:,:,1,:),[1,2,3],'omitnan')),squeeze(std(noise_masker_plot_data_pfc(:,:,1,:),[],[1,2,3],'omitnan'))./sqrt(num_subjects-1), 'lineProps',colors(1) )
shadedErrorBar(time,squeeze(mean(noise_masker_plot_data_pfc(:,:,2,:),[1,2,3],'omitnan')),squeeze(std(noise_masker_plot_data_pfc(:,:,2,:),[],[1,2,3],'omitnan'))./sqrt(num_subjects-1), 'lineProps',colors(2)  )
ylim([ymin,ymax])
xlim([epoch_start_time,epoch_end_time])
set(gca,'fontsize',14)


nexttile
hold on
shadedErrorBar(time,squeeze(mean(noise_masker_plot_data_pfc(:,:,3,:),[1,2,3],'omitnan')),squeeze(std(noise_masker_plot_data_pfc(:,:,3,:),[],[1,2,3],'omitnan'))./sqrt(num_subjects-1), 'lineProps',colors(3)  )
shadedErrorBar(time,squeeze(mean(noise_masker_plot_data_pfc(:,:,4,:),[1,2,3],'omitnan')),squeeze(std(noise_masker_plot_data_pfc(:,:,4,:),[],[1,2,3],'omitnan'))./sqrt(num_subjects-1), 'lineProps',colors(4)  )
ylim([ymin,ymax])
xlim([epoch_start_time,epoch_end_time])
set(gca,'fontsize',14)


% STG, noise Masker
nexttile
hold on
shadedErrorBar(time,squeeze(mean(noise_masker_plot_data_stg(:,:,1,:),[1,2,3],'omitnan')),squeeze(std(noise_masker_plot_data_stg(:,:,1,:),[],[1,2,3],'omitnan'))./sqrt(num_subjects-1), 'lineProps',colors(1) )
shadedErrorBar(time,squeeze(mean(noise_masker_plot_data_stg(:,:,2,:),[1,2,3],'omitnan')),squeeze(std(noise_masker_plot_data_stg(:,:,2,:),[],[1,2,3],'omitnan'))./sqrt(num_subjects-1), 'lineProps',colors(2)  )
ylim([ymin,ymax])
xlim([epoch_start_time,epoch_end_time])
set(gca,'fontsize',14)

nexttile
hold on
shadedErrorBar(time,squeeze(mean(noise_masker_plot_data_stg(:,:,3,:),[1,2,3],'omitnan')),squeeze(std(noise_masker_plot_data_stg(:,:,3,:),[],[1,2,3],'omitnan'))./sqrt(num_subjects-1), 'lineProps',colors(3)  )
shadedErrorBar(time,squeeze(mean(noise_masker_plot_data_stg(:,:,4,:),[1,2,3],'omitnan')),squeeze(std(noise_masker_plot_data_stg(:,:,4,:),[],[1,2,3],'omitnan'))./sqrt(num_subjects-1), 'lineProps',colors(4)  )
ylim([ymin,ymax])
xlim([epoch_start_time,epoch_end_time])
set(gca,'fontsize',14)
