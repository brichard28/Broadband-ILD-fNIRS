%% srm_nirs_eeg_plot_mean_hbo

% load in data
speech_masker_data = readtable("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_speech_masker.csv",'format','auto');
noise_masker_data = readtable("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_noise_masker.csv",'format','auto');
speech_masker_data = speech_masker_data(2:end,:);
noise_masker_data = noise_masker_data(2:end,:);

speech_masker_data.Properties.VariableNames  = {'S','Channel','ITD50', 'ITD500', 'ILD70n', 'ILD10'};
noise_masker_data.Properties.VariableNames  = {'S','Channel','ITD50', 'ITD500', 'ILD70n', 'ILD10'};

speech_masker_data = stack(speech_masker_data,{["ITD50" "ITD500" "ILD70n" "ILD10"]}, ...
    "NewDataVariableName",["MeanHbO"], "IndexVariableName","Spatialization");
noise_masker_data = stack(noise_masker_data,{["ITD50" "ITD500" "ILD70n" "ILD10"]}, ...
    "NewDataVariableName",["MeanHbO"], "IndexVariableName","Spatialization");

% concatenate speech and noise
all_data = [speech_masker_data;noise_masker_data];
all_data.Masker = [repmat("speech",height(speech_masker_data),1);repmat("noise",height(noise_masker_data),1)];
for irow = 1:height(all_data)
    this_channel = all_data(irow,'Channel');
    if ismember(str2double(this_channel{1,1}),[0,1,2,3,4,5])
        all_data(irow,'Roi') = {'pfc'};
    else
        all_data(irow,'Roi') = {'stg'};
    end
end

subject_ID = unique(all_data.S);
%% Compare across masker in each channel
figure;
tiledlayout(2,7)
for ichannel = 0:13
    nexttile
    hold on
    speech_masker_plot_data = all_data.MeanHbO(all_data.Masker == "speech" & str2double(all_data.Channel) == ichannel);
    noise_masker_plot_data = all_data.MeanHbO(all_data.Masker == "noise" & str2double(all_data.Channel) == ichannel);
    line_color_rgb = 0.8;

    plot([1:2],[speech_masker_plot_data,noise_masker_plot_data], 'Color', [line_color_rgb line_color_rgb line_color_rgb])
    espeech = errorbar(1,mean(speech_masker_plot_data,1,'omitnan'),std(speech_masker_plot_data,[],1,'omitnan')./(sqrt(size(subject_ID,1) -1)));
    enoise = errorbar(2,mean(noise_masker_plot_data,1,'omitnan'),std(noise_masker_plot_data,[],1,'omitnan')./(sqrt(size(subject_ID,1) -1)));

    espeech.Marker = 'o';
    espeeech.MarkerSize = 4;
    espeech.MarkerFaceColor = 'm';
    espeech.Color = 'm';
    espeech.CapSize = 2;
    espeech.LineWidth = 2;
    enoise.Marker = 'o';
    enoise.MarkerSize = 4;
    enoise.MarkerFaceColor = 'c';
    enoise.Color = 'c';
    enoise.CapSize = 2;
    enoise.LineWidth = 2;

    xlim([0.75, 2.25])
    xticks([1,2])
end

xlabel("Masker",'FontSize',18,'FontWeight','bold')
xticklabels({"Speech","Noise"})
ylabel('Mean \DeltaHbO (\muM)','FontSize',18,'FontWeight','bold')
%% Compare across ROI

% stg_plot_data = all_data.MeanHbO(all_data.Roi == "stg");
% pfc_plot_data = all_data.MeanHbO(all_data.Roi == "pfc");
% line_color_rgb = 0.8;
% figure;
% hold on
% plot([1:2],[pfc_plot_data, stg_plot_data], 'Color', [line_color_rgb line_color_rgb line_color_rgb])
% epfc = errorbar(2,mean(pfc_plot_data,'omitnan'),std(pfc_plot_data,[],1,'omitnan')./(sqrt(size(subject_ID,1) -1)));
% estg = errorbar(1,mean(stg_plot_data,1,'omitnan'),std(stg_plot_data,[],1,'omitnan')./(sqrt(size(subject_ID,1) -1)));
%
%
% xlim([0.75, 2.25])
% xticks([1,2])
% xlabel("ROI",'FontSize',18,'FontWeight','bold')
% xticklabels({"PFC","STG"})
% ylabel('Mean \DeltaHbO (\muM)','FontSize',18,'FontWeight','bold')

%% Compare across spatialization and masker (2 panel, each ROI)
% PFC

for ichannel = 0:13
figure;
tiledlayout(1,2,'TileSpacing','compact');
nexttile
hold on
% Across spatialization, speech masker
this_data = all_data(all_data.Masker == "speech" & str2double(all_data.Channel) == ichannel, :);
plot([1:2],[this_data.MeanHbO(this_data.Spatialization == "ITD50"),this_data.MeanHbO(this_data.Spatialization == "ITD500")],'Color',[line_color_rgb line_color_rgb line_color_rgb])
plot([2.5:3.5],[this_data.MeanHbO(this_data.Spatialization == "ILD70n"),this_data.MeanHbO(this_data.Spatialization == "ILD10")],'Color',[line_color_rgb line_color_rgb line_color_rgb])

e1 = errorbar(1,mean(this_data.MeanHbO(this_data.Spatialization == "ITD50"),1,'omitnan'),std(this_data.MeanHbO(this_data.Spatialization == "ITD50"),[],1,'omitnan')./(sqrt(size(subject_ID,1))-1));
e2 = errorbar(2,mean(this_data.MeanHbO(this_data.Spatialization == "ITD500"),1,'omitnan'),std(this_data.MeanHbO(this_data.Spatialization == "ITD500"),[],1,'omitnan')./(sqrt(size(subject_ID,1))-1));
e3 = errorbar(2.5,mean(this_data.MeanHbO(this_data.Spatialization == "ILD70n"),1,'omitnan'),std(this_data.MeanHbO(this_data.Spatialization == "ILD70n"),[],1,'omitnan')./(sqrt(size(subject_ID,1))-1));
e4 = errorbar(3.5,mean(this_data.MeanHbO(this_data.Spatialization == "ILD10"),1,'omitnan'),std(this_data.MeanHbO(this_data.Spatialization == "ILD10"),[],1,'omitnan')./(sqrt(size(subject_ID,1))-1));

e1.Marker = 'o';
e1.MarkerSize = 4;
e1.MarkerFaceColor = 'r';
e1.Color = 'r';
e1.CapSize = 2;
e1.LineWidth = 2;
e2.Marker = 'o';
e2.MarkerSize = 4;
e2.MarkerFaceColor = 'blue';
e2.Color = 'blue';
e2.CapSize = 2;
e2.LineWidth = 2;
e3.Marker = 'o';
e3.MarkerSize = 4;
e3.MarkerFaceColor = 'm';
e3.Color = 'm';
e3.CapSize = 2;
e3.LineWidth = 2;
e4.Marker = 'o';
e4.MarkerSize = 4;
e4.MarkerFaceColor = 'g';
e4.Color = 'g';
e4.CapSize = 2;
e4.LineWidth = 2;
ax = gca;
ax.LineWidth = 1;
ax.XTickLabel = cellfun(@(a) ['\bf{' a '}'], ax.XTickLabel, 'UniformOutput',false);

xticks([1,2,2.5,3.5])
xlim([0.75 3.75])
ylim([-0.1, 0.4])
xticklabels({'Small ITD','Large ITD','Natural ILD','Broadband ILD'})
ylabel('Mean \DeltaHbO (\muM)','FontSize',18,'FontWeight','bold')
title('Speech Masker','FontSize',18)


nexttile
hold on
% Across spatialization, speech masker
this_data = all_data(all_data.Masker == "noise" & str2double(all_data.Channel) == ichannel, :);
plot([1:2],[this_data.MeanHbO(this_data.Spatialization == "ITD50"),this_data.MeanHbO(this_data.Spatialization == "ITD500")],'Color',[line_color_rgb line_color_rgb line_color_rgb])
plot([2.5:3.5],[this_data.MeanHbO(this_data.Spatialization == "ILD70n"),this_data.MeanHbO(this_data.Spatialization == "ILD10")],'Color',[line_color_rgb line_color_rgb line_color_rgb])

e1 = errorbar(1,mean(this_data.MeanHbO(this_data.Spatialization == "ITD50"),1,'omitnan'),std(this_data.MeanHbO(this_data.Spatialization == "ITD50"),[],1,'omitnan')./(sqrt(size(subject_ID,1))-1));
e2 = errorbar(2,mean(this_data.MeanHbO(this_data.Spatialization == "ITD500"),1,'omitnan'),std(this_data.MeanHbO(this_data.Spatialization == "ITD500"),[],1,'omitnan')./(sqrt(size(subject_ID,1))-1));
e3 = errorbar(2.5,mean(this_data.MeanHbO(this_data.Spatialization == "ILD70n"),1,'omitnan'),std(this_data.MeanHbO(this_data.Spatialization == "ILD70n"),[],1,'omitnan')./(sqrt(size(subject_ID,1))-1));
e4 = errorbar(3.5,mean(this_data.MeanHbO(this_data.Spatialization == "ILD10"),1,'omitnan'),std(this_data.MeanHbO(this_data.Spatialization == "ILD10"),[],1,'omitnan')./(sqrt(size(subject_ID,1))-1));

e1.Marker = 'o';
e1.MarkerSize = 4;
e1.MarkerFaceColor = 'r';
e1.Color = 'r';
e1.CapSize = 2;
e1.LineWidth = 2;
e2.Marker = 'o';
e2.MarkerSize = 4;
e2.MarkerFaceColor = 'blue';
e2.Color = 'blue';
e2.CapSize = 2;
e2.LineWidth = 2;
e3.Marker = 'o';
e3.MarkerSize = 4;
e3.MarkerFaceColor = 'm';
e3.Color = 'm';
e3.CapSize = 2;
e3.LineWidth = 2;
e4.Marker = 'o';
e4.MarkerSize = 4;
e4.MarkerFaceColor = 'g';
e4.Color = 'g';
e4.CapSize = 2;
e4.LineWidth = 2;
ax = gca;
ax.LineWidth = 1;
ax.XTickLabel = cellfun(@(a) ['\bf{' a '}'], ax.XTickLabel, 'UniformOutput',false);

xticks([1,2,2.5,3.5])
xlim([0.75 3.75])
ylim([-0.1, 0.4])
xticklabels({'Small ITD','Large ITD','Natural ILD','Broadband ILD'})
ylabel('Mean \DeltaHbO (\muM)','FontSize',18,'FontWeight','bold')
title('Noise Masker','FontSize',18)


sgtitle(ichannel)
end