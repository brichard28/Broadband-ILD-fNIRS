%% BNR_make_topoplots

% Make topographic plots for focal versus broad

subject_tags = {'018XC','019XX','021XC','022XX','023XX','025MX','026MC','027MC',...
    '028MC','029MX','030MX','031MC',...
    '031XX','032MX','034MX','035MX','039XC','040MC',...
    '051XC','053XX','055XX','061XX','063XX','064XC','067XX',...
    '068XX','070XX','073XX','075XC','076XX','080XX','084MX',...
    '091MC','093XX','096MC','098MX',...
    '100MX','101XC','102XX',...
    '103MX','105MX','110XX',...
    '112MX','200XC'};
is_adhd = [0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,1,1,1,1,1,1,0,1,1,1,1,1,0];

load('Alpha_Powers_BNR_ 2022_5_5_7_3.mat','all_alpha_powers_baselined','times','freqs')

for isubject = 1:length(subject_tags) - 1
    curr_subj_tag = split(string(subject_tags(isubject)),'');
    curr_subj_tag = curr_subj_tag(2:6)';
    load(join(string(['AUD_ASA',curr_subj_tag(1:3),'_',curr_subj_tag(4:5),'_1-50_DATA.mat']),''),'SCORE');
    correct_trials = find(SCORE.hits(:)'==1);
    delay_1s = SCORE.delay == 1;
    F_N_trials = find((SCORE.F_N(:)'.*SCORE.hits(:)')==1);
    F_S1_trials = find((SCORE.F_S1(:)'.*SCORE.hits(:)')==1);
    F_S2_trials = find((SCORE.F_S2(:)'.*SCORE.hits(:)')==1);
    B_N_trials = find((SCORE.B_N(:)'.*SCORE.hits(:)')==1);
    B_S1_trials = find((SCORE.B_S1(:)'.*SCORE.hits(:)')==1);
    B_S2_trials = find((SCORE.B_S2(:)'.*SCORE.hits(:)')==1);
    FOCAL_trials = [F_N_trials,F_S1_trials,F_S2_trials];
    BROAD_trials = [B_N_trials,B_S1_trials,B_S2_trials];

    % all focal
    current_FOCAL_power = all_alpha_powers_baselined(isubject,FOCAL_trials,:,:);
    FOCAL_alpha_allsubjects(isubject,:,:) = squeeze(nanmean(current_FOCAL_power,2)); % mean across focal trials

    % F_N
    current_F_N_power = all_alpha_powers_baselined(isubject,F_N_trials,:,:);
    F_N_alpha_allsubjects(isubject,:,:) = squeeze(nanmean(current_F_N_power,2)); % mean across focal trials

    % F_S1
    current_F_S1_power = all_alpha_powers_baselined(isubject,F_S1_trials,:,:);
    F_S1_alpha_allsubjects(isubject,:,:) = squeeze(nanmean(current_F_S1_power,2)); % mean across focal trials

    % F_S2
    current_F_S2_power = all_alpha_powers_baselined(isubject,F_S2_trials,:,:);
    F_S2_alpha_allsubjects(isubject,:,:) = squeeze(nanmean(current_F_S2_power,2)); % mean across focal trials

   
    % all broad
    current_BROAD_power = all_alpha_powers_baselined(isubject,BROAD_trials,:,:);
    BROAD_alpha_allsubjects(isubject,:,:) = squeeze(nanmean(current_BROAD_power,2)); % mean across broad trials

    % B_N 
    current_B_N_power = all_alpha_powers_baselined(isubject,B_N_trials,:,:);
    B_N_alpha_allsubjects(isubject,:,:) = squeeze(nanmean(current_B_N_power,2)); % mean across focal trials

    % B_S1 
    current_B_S1_power = all_alpha_powers_baselined(isubject,B_S1_trials,:,:);
    B_S1_alpha_allsubjects(isubject,:,:) = squeeze(nanmean(current_B_S1_power,2)); % mean across focal trials

    % B_S2 
    current_B_S2_power = all_alpha_powers_baselined(isubject,B_S2_trials,:,:);
    B_S2_alpha_allsubjects(isubject,:,:) = squeeze(nanmean(current_B_S2_power,2)); % mean across focal trials

end


%% Make Focal topoplot
figure;
FOCAL_alpha = squeeze(mean(FOCAL_alpha_allsubjects,[1,3])); % average across subject and time
topoplot(FOCAL_alpha,'chanlocs_64_chan.locs');
colorbar
%caxis([-1 1])
title('Focal Only Average Alpha Power')

%% Make F_N topoplot
figure;
F_N_alpha = squeeze(nanmean(F_N_alpha_allsubjects,[1,3])); % average across subject and time
topoplot(F_N_alpha,'chanlocs_64_chan.locs');
colorbar
title('F_N Average Alpha Power')

%% Make F_S1 topoplot
figure;
F_S1_alpha = squeeze(nanmean(F_S1_alpha_allsubjects,[1,3])); % average across subject and time
topoplot(F_S1_alpha,'chanlocs_64_chan.locs');
colorbar
%caxis([-1 1])
title('F_S1 Average Alpha Power')

%% Make F_S2 topoplot
figure;
F_S2_alpha = squeeze(nanmean(F_S2_alpha_allsubjects,[1,3])); % average across subject and time
topoplot(F_S2_alpha,'chanlocs_64_chan.locs');
colorbar
%caxis([-1 1])

title('F_S2 Average Alpha Power')


%% Make Broad topoplot
figure;
BROAD_alpha = squeeze(nanmean(BROAD_alpha_allsubjects,[1,3])); % average across subject and time
topoplot(BROAD_alpha,'chanlocs_64_chan.locs');
colorbar
%caxis([-1 1])

title('Broad Only Average Alpha Power')

%% Make B_N topoplot
figure;
B_N_alpha = squeeze(nanmean(B_N_alpha_allsubjects,[1,3])); % average across subject and time
topoplot(B_N_alpha,'chanlocs_64_chan.locs');
colorbar
%caxis([-1 1])

title('B_N Average Alpha Power')


%% Make B_S1 topoplot
figure;
B_S1_alpha = squeeze(nanmean(B_S1_alpha_allsubjects,[1,3])); % average across subject and time
topoplot(B_S1_alpha,'chanlocs_64_chan.locs');
colorbar
%caxis([-1 1])

title('B_S1 Average Alpha Power')

%% Make B_S2 topoplot
figure;
B_S2_alpha = squeeze(nanmean(B_S2_alpha_allsubjects,[1,3])); % average across subject and time
topoplot(B_S2_alpha,'chanlocs_64_chan.locs');
colorbar
%caxis([-1 1])

title('B_S2 Average Alpha Power')


%% Make Focal minus Broad topoplot
figure;
DIFF_alpha = squeeze(nanmean(FOCAL_alpha_allsubjects,[1,3])) - squeeze(nanmean(BROAD_alpha_allsubjects,[1,3])); % average across subject and time
topoplot(DIFF_alpha,'chanlocs_64_chan.locs');
colorbar
%caxis([-1 1])

title('Focal Minus Broad Average Alpha Power')

%% Make F_N minus B_N topoplot
figure;
DIFF_alpha_N = squeeze(nanmean(F_N_alpha_allsubjects,[1,3])) - squeeze(nanmean(B_N_alpha_allsubjects,[1,3])); % average across subject and time
topoplot(DIFF_alpha_N,'chanlocs_64_chan.locs');
colorbar
%caxis([-1 1])

title('F_N Minus B_N Average Alpha Power','Interpreter','None')

%% Make F_S1 minus B_S1 topoplot
figure;
DIFF_alpha_S1 = squeeze(nanmean(F_S1_alpha_allsubjects,[1,3])) - squeeze(nanmean(B_S1_alpha_allsubjects,[1,3])); % average across subject and time
topoplot(DIFF_alpha_S1,'chanlocs_64_chan.locs');
colorbar
%caxis([-1 1])

title('F_S1 Minus B_S1 Average Alpha Power','Interpreter','None')

%% Make F_S2 minus B_S2 topoplot
figure;
DIFF_alpha_S2 = squeeze(nanmean(F_S2_alpha_allsubjects,[1,3])) - squeeze(nanmean(B_S2_alpha_allsubjects,[1,3])); % average across subject and time
topoplot(DIFF_alpha_S2,'chanlocs_64_chan.locs');
colorbar
%caxis([-1 1])

title('F_S2 Minus B_S2 Average Alpha Power','Interpreter','None')

%% Movie sample rate
srate = length(times)/7;
% %% Make a Focal movie
% figure;
% data = squeeze(mean(FOCAL_alpha_allsubjects,1));
% [FOCAL_Movie,FOCAL_Colormap] = eegmovie(data,srate,'chanlocs_64_chan.locs'); %% Is matlab forcing iso-lines to be same values?? Double check
% v = VideoWriter('All_Focal_Movie','Archival');
% open(v)
% writeVideo(v,FOCAL_Movie)
% 
% %% Make a Broad movie
% figure;
% data = squeeze(mean(BROAD_alpha_allsubjects,1));
% [BROAD_Movie,BROAD_Colormap] = eegmovie(data,srate,'chanlocs_64_chan.locs');
% v = VideoWriter('All_Broad_Movie','Archival');
% open(v)
% writeVideo(v, BROAD_Movie)
% 
% %% F_N movie 
% figure;
% data = squeeze(mean(F_N_alpha_allsubjects,1));
% [F_N_Movie,F_N_Colormap] = eegmovie(data,srate,'chanlocs_64_chan.locs');
% v = VideoWriter('F_N_Movie','Archival');
% open(v)
% writeVideo(v,F_N_Movie)
% 
% 
% %% F_S1 movie
% figure;
% data = squeeze(mean(F_S1_alpha_allsubjects,1));
% [F_S1_Movie,F_S1_Colormap] = eegmovie(data,srate,'chanlocs_64_chan.locs');
% v = VideoWriter('F_S1_Movie','Archival');
% open(v)
% writeVideo(v,F_S1_Movie)
% 
% 
% %% F_S2 movie
% figure;
% data = squeeze(mean(F_S2_alpha_allsubjects,1));
% [F_S2_Movie,F_S2_Colormap] = eegmovie(data,srate,'chanlocs_64_chan.locs');
% v = VideoWriter('F_S2_Movie','Archival');
% open(v)
% writeVideo(v,F_S2_Movie)
% 
% %% B_N movie 
% figure;
% data = squeeze(mean(B_N_alpha_allsubjects,1));
% [B_N_Movie,B_N_Colormap] = eegmovie(data,srate,'chanlocs_64_chan.locs');
% v = VideoWriter('B_N_Movie','Archival');
% open(v)
% writeVideo(v,B_N_Movie)
% 
% %% B_S1 movie
% figure;
% data = squeeze(mean(B_S1_alpha_allsubjects,1));
% [B_S1_Movie,B_S1_Colormap] = eegmovie(data,srate,'chanlocs_64_chan.locs');
% v = VideoWriter('B_S1_Movie','Archival');
% open(v)
% writeVideo(v,B_S1_Movie)
% 
% 
% %% B_S2 movie
% figure;
% data = squeeze(mean(B_S2_alpha_allsubjects,1));
% [B_S2_Movie,B_S2_Colormap] = eegmovie(data,srate,'chanlocs_64_chan.locs');
% v = VideoWriter('B_S2_Movie','Archival');
% open(v)
% writeVideo(v,B_S2_Movie)

% %% Difference movie
% figure;
% DIFF_alpha_allsubjects = FOCAL_alpha_allsubjects - BROAD_alpha_allsubjects;
% data = squeeze(mean(DIFF_alpha_allsubjects,1));
% [Difference_movie,Difference_Colormap] = eegmovie(data,srate,'chanlocs_64_chan.locs');
% colorbar;
%  v = VideoWriter('Difference Movie','Archival');
%  open(v)
%  writeVideo(v,Difference_movie)
% close(v)