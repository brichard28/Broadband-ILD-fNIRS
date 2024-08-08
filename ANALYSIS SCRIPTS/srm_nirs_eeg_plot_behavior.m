%% srm_nirs_eeg_plot_behavior

load('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_Behavior_Results.mat')
subject_ID = 1:size(all_hit_rates_collapsed,2);
subject_ID = subject_ID';

% NOTES:
% change from red and blue between condition pairs --> thats for fNIRS

% left column = speech masker, right column = noise masker
% top row = color word hit rate, bottom row = object rate
figure;
tiledlayout(2, 3, 'TileSpacing','tight');
% TOP ROW: Speech Masker Hit, Speech Masker Object, Speech Masker FA,
% Speech Masker D- Prime

%% Speech Masker Hit Rate
nexttile
hold on
for i = 1:size(all_hit_rates_collapsed,2)
    if all_hit_rates_collapsed(5,i) > all_hit_rates_collapsed(6,i) % plot in red
        plot([1:2],all_hit_rates_collapsed(5:6,i),'Color',[0.8 0.8 0.8])
    else % plot in blue
        plot([1:2],all_hit_rates_collapsed(5:6,i),'Color',[0.8 0.8 0.8])
    end
    if all_hit_rates_collapsed(7,i) > all_hit_rates_collapsed(8,i) % plot in red
       plot([2.5:3.5],all_hit_rates_collapsed(7:8,i),'Color',[0.8 0.8 0.8])
    else % plot in blue
        plot([2.5:3.5],all_hit_rates_collapsed(7:8,i),'Color',[0.8 0.8 0.8])
    end
end
e1 = errorbar(1,mean(all_hit_rates_collapsed(5,:),2,'omitnan'),std(all_hit_rates_collapsed(5,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e2 = errorbar(2,mean(all_hit_rates_collapsed(6,:),2,'omitnan'),std(all_hit_rates_collapsed(6,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e3 = errorbar(2.5,mean(all_hit_rates_collapsed(7,:),2,'omitnan'),std(all_hit_rates_collapsed(7,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e4 = errorbar(3.5,mean(all_hit_rates_collapsed(8,:),2,'omitnan'),std(all_hit_rates_collapsed(8,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));

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
ylabel('Speech \newline Masker','FontSize',18,'FontWeight','bold','Rotation',0)
ylim([0,1])
title('Hit Rate (%)','FontSize',18)
set(gca,'XTickLabel',[]);

% significance bars
line([1,2],[0.7,0.7],'Color','k','LineWidth',1.5)
text(1.5,0.72,'***','FontSize',18)
line([1,2.5],[0.8,0.8],'Color','k','LineWidth',1.5)
text(1.75,0.82,'***','FontSize',18)
line([1,3.5],[0.9,0.9],'Color','k','LineWidth',1.5)
text(2.25,0.92,'***','FontSize',18)


%% Speech Masker Object Rate
% nexttile
% hold on
% for i = 1:size(all_object_rates_collapsed,2)
%     if all_object_rates_collapsed(5,i) > all_object_rates_collapsed(6,i) % plot in red
%         plot([1:2],all_object_rates_collapsed(5:6,i),'Color',[0.8 0.8 0.8])
%     else % plot in blue
%         plot([1:2],all_object_rates_collapsed(5:6,i),'Color',[0.8 0.8 0.8])
%     end
%     if all_object_rates_collapsed(7,i) > all_object_rates_collapsed(8,i) % plot in red
%        plot([2.5:3.5],all_object_rates_collapsed(7:8,i),'Color',[0.8 0.8 0.8])
%     else % plot in blue
%         plot([2.5:3.5],all_object_rates_collapsed(7:8,i),'Color',[0.8 0.8 0.8])
%     end
% end
% e1 = errorbar(1,mean(all_object_rates_collapsed(5,:),2,'omitnan'),std(all_object_rates_collapsed(5,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
% e2 = errorbar(2,mean(all_object_rates_collapsed(6,:),2,'omitnan'),std(all_object_rates_collapsed(6,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
% e3 = errorbar(2.5,mean(all_object_rates_collapsed(7,:),2,'omitnan'),std(all_object_rates_collapsed(7,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
% e4 = errorbar(3.5,mean(all_object_rates_collapsed(8,:),2,'omitnan'),std(all_object_rates_collapsed(8,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
% 
% e1.Marker = 'o';
% e1.MarkerSize = 4;
% e1.MarkerFaceColor = 'r';
% e1.Color = 'r';
% e1.CapSize = 2;
% e1.LineWidth = 2;
% e2.Marker = 'o';
% e2.MarkerSize = 4;
% e2.MarkerFaceColor = 'blue';
% e2.Color = 'blue';
% e2.CapSize = 2;
% e2.LineWidth = 2;
% e3.Marker = 'o';
% e3.MarkerSize = 4;
% e3.MarkerFaceColor = 'm';
% e3.Color = 'm';
% e3.CapSize = 2;
% e3.LineWidth = 2;
% e4.Marker = 'o';
% e4.MarkerSize = 4;
% e4.MarkerFaceColor = 'g';
% e4.Color = 'g';
% e4.CapSize = 2;
% e4.LineWidth = 2;
% ax = gca;
% ax.LineWidth = 1;
% ax.XTickLabel = cellfun(@(a) ['\bf{' a '}'], ax.XTickLabel, 'UniformOutput',false);
% 
% title('Object Rate (%)','FontSize',18,'FontWeight','bold')
% xticks([1,2,2.5,3.5])
% xlim([0.75 3.75])
% ylim([0,1])
% 
% % significance bars
% line([1,2.5],[0.1,0.1],'Color','k','LineWidth',1.5)
% text(1.65,0.11,'***','FontSize',14)
% line([2,2.5],[0.125,0.125],'Color','k','LineWidth',1.5)
% text(2.15,0.135,'***','FontSize',14)
% line([2.5,3.5],[0.15,0.15],'Color','k','LineWidth',1.5)
% text(2.9,0.16,'***','FontSize',14)

%% Speech Masker FA Rate
nexttile
hold on
for i = 1:size(all_FA_rates_collapsed,2)
    if all_FA_rates_collapsed(5,i) > all_FA_rates_collapsed(6,i) % plot in red
        plot([1:2],all_FA_rates_collapsed(5:6,i),'Color',[0.8 0.8 0.8])
    else % plot in blue
        plot([1:2],all_FA_rates_collapsed(5:6,i),'Color',[0.8 0.8 0.8])
    end
    if all_FA_rates_collapsed(7,i) > all_FA_rates_collapsed(8,i) % plot in red
       plot([2.5:3.5],all_FA_rates_collapsed(7:8,i),'Color',[0.8 0.8 0.8])
    else % plot in blue
        plot([2.5:3.5],all_FA_rates_collapsed(7:8,i),'Color',[0.8 0.8 0.8])
    end
end
e1 = errorbar(1,mean(all_FA_rates_collapsed(5,:),2,'omitnan'),std(all_FA_rates_collapsed(5,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e2 = errorbar(2,mean(all_FA_rates_collapsed(6,:),2,'omitnan'),std(all_FA_rates_collapsed(6,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e3 = errorbar(2.5,mean(all_FA_rates_collapsed(7,:),2,'omitnan'),std(all_FA_rates_collapsed(7,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e4 = errorbar(3.5,mean(all_FA_rates_collapsed(8,:),2,'omitnan'),std(all_FA_rates_collapsed(8,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));

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
xticklabels({'Small ITD','Large ITD','Natural ILD','Broadband ILD'})
title('FA Rate (%)','FontSize',18,'FontWeight','bold')

% significance bars
line([1,2.5],[0.4,0.4],'Color','k','LineWidth',1.5)
text(1.65,0.42,'***','FontSize',18)
line([2,2.5],[0.5,0.5],'Color','k','LineWidth',1.5)
text(2.15,0.52,'***','FontSize',18)
line([2.5,3.5],[0.6,0.6],'Color','k','LineWidth',1.5)
text(2.9,0.62,'***','FontSize',18)
ylim([0,1])
set(gca,'YTickLabel',[]);

%% Speech Masker d'
nexttile
hold on
for i = 1:size(d_primes_collapsed,2)
    if d_primes_collapsed(5,i) > d_primes_collapsed(6,i) % plot in red
        plot([1:2],d_primes_collapsed(5:6,i),'Color',[0.8 0.8 0.8])
    else % plot in blue
        plot([1:2],d_primes_collapsed(5:6,i),'Color',[0.8 0.8 0.8])
    end
    if d_primes_collapsed(7,i) > d_primes_collapsed(8,i) % plot in red
       plot([2.5:3.5],d_primes_collapsed(7:8,i),'Color',[0.8 0.8 0.8])
    else % plot in blue
        plot([2.5:3.5],d_primes_collapsed(7:8,i),'Color',[0.8 0.8 0.8])
    end
end
e1 = errorbar(1,mean(d_primes_collapsed(5,:),2,'omitnan'),std(d_primes_collapsed(5,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e2 = errorbar(2,mean(d_primes_collapsed(6,:),2,'omitnan'),std(d_primes_collapsed(6,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e3 = errorbar(2.5,mean(d_primes_collapsed(7,:),2,'omitnan'),std(d_primes_collapsed(7,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e4 = errorbar(3.5,mean(d_primes_collapsed(8,:),2,'omitnan'),std(d_primes_collapsed(8,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));

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
xticklabels({'Small ITD','Large ITD','Natural ILD','Broadband ILD'})
title("Behavioral\newlineSensitivity (d')",'FontSize',18,'FontWeight','bold')


%% Noise Masker Hit Rate
nexttile
hold on
for i = 1:size(all_hit_rates_collapsed,2)
    if all_hit_rates_collapsed(1,i) > all_hit_rates_collapsed(2,i) % plot in red
        plot([1:2],all_hit_rates_collapsed(1:2,i),'Color',[0.8 0.8 0.8])
    else % plot in blue
        plot([1:2],all_hit_rates_collapsed(1:2,i),'Color',[0.8 0.8 0.8])
    end
    if all_hit_rates_collapsed(3,i) > all_hit_rates_collapsed(4,i) % plot in red
       plot([2.5:3.5],all_hit_rates_collapsed(3:4,i),'Color',[0.8 0.8 0.8])
    else % plot in blue
        plot([2.5:3.5],all_hit_rates_collapsed(3:4,i),'Color',[0.8 0.8 0.8])
    end
end
e1 = errorbar(1,mean(all_hit_rates_collapsed(1,:),2,'omitnan'),std(all_hit_rates_collapsed(1,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e2 = errorbar(2,mean(all_hit_rates_collapsed(2,:),2,'omitnan'),std(all_hit_rates_collapsed(2,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e3 = errorbar(2.5,mean(all_hit_rates_collapsed(3,:),2,'omitnan'),std(all_hit_rates_collapsed(3,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
e4 = errorbar(3.5,mean(all_hit_rates_collapsed(4,:),2,'omitnan'),std(all_hit_rates_collapsed(4,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));

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
xticklabels({'Small ITD','Large ITD','Natural ILD','Broadband ILD'})
xlim([0.75 3.75])
ylabel('Noise\newlineMasker','FontSize',18,'FontWeight','bold','Rotation',0)
ylim([0,1])

%% Noise Masker Object Rate
% nexttile
% hold on
% for i = 1:size(all_object_rates_collapsed,2)
%     if all_object_rates_collapsed(1,i) > all_object_rates_collapsed(2,i) % plot in red
%         plot([1:2],all_object_rates_collapsed(1:2,i),'Color',[0.8 0.8 0.8])
%     else % plot in blue
%         plot([1:2],all_object_rates_collapsed(1:2,i),'Color',[0.8 0.8 0.8])
%     end
%     if all_object_rates_collapsed(3,i) > all_object_rates_collapsed(4,i) % plot in red
%        plot([2.5:3.5],all_object_rates_collapsed(3:4,i),'Color',[0.8 0.8 0.8])
%     else % plot in blue
%         plot([2.5:3.5],all_object_rates_collapsed(3:4,i),'Color',[0.8 0.8 0.8])
%     end
% end
% e1 = errorbar(1,mean(all_object_rates_collapsed(1,:),2,'omitnan'),std(all_object_rates_collapsed(1,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
% e2 = errorbar(2,mean(all_object_rates_collapsed(2,:),2,'omitnan'),std(all_object_rates_collapsed(2,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
% e3 = errorbar(2.5,mean(all_object_rates_collapsed(3,:),2,'omitnan'),std(all_object_rates_collapsed(3,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
% e4 = errorbar(3.5,mean(all_object_rates_collapsed(4,:),2,'omitnan'),std(all_object_rates_collapsed(4,:),[],2,'omitnan')./(sqrt(size(subject_ID,1))-1));
% 
% e1.Marker = 'o';
% e1.MarkerSize = 4;
% e1.MarkerFaceColor = 'r';
% e1.Color = 'r';
% e1.CapSize = 2;
% e1.LineWidth = 2;
% e2.Marker = 'o';
% e2.MarkerSize = 4;
% e2.MarkerFaceColor = 'blue';
% e2.Color = 'blue';
% e2.CapSize = 2;
% e2.LineWidth = 2;
% e3.Marker = 'o';
% e3.MarkerSize = 4;
% e3.MarkerFaceColor = 'm';
% e3.Color = 'm';
% e3.CapSize = 2;
% e3.LineWidth = 2;
% e4.Marker = 'o';
% e4.MarkerSize = 4;
% e4.MarkerFaceColor = 'g';
% e4.Color = 'g';
% e4.CapSize = 2;
% e4.LineWidth = 2;
% ax = gca;
% ax.LineWidth = 1;
% ax.XTickLabel = cellfun(@(a) ['\bf{' a '}'], ax.XTickLabel, 'UniformOutput',false);
% 
% xticks([1,2,2.5,3.5])
% xlim([0.75 3.75])
% xticklabels({'Small ITD','Large ITD','Natural ILD','Broadband ILD'})
% ylim([0,1])







%% Plot hit-rate vs d-prime for speech masker
% figure;
% tiledlayout(2,2)
% hold on
% colors = {'r','b','m','g'};
% leg_labels = {'Small ITD','Large ITD','Natural ILD','Broadband ILD'};
% 
% for icondition = 1:4
%     nexttile
%     hold on
%     curr_d_primes = d_primes_collapsed(icondition+4,:);
%     curr_hit_rates = all_hit_rates_collapsed(icondition+4,:);
%     scatter(curr_d_primes,curr_hit_rates,'filled',string(colors(icondition)))
%     [R,P] = corrcoef(curr_d_primes,curr_hit_rates);
%     [p_fit,S] = polyfit(curr_d_primes,curr_hit_rates,1);
%     [fit_line,delta] = polyval(p_fit,curr_d_primes,S);
%     R_squared = 1 - (S.normr/norm(curr_hit_rates - mean(curr_hit_rates)))^2;
%     plot(curr_d_primes,fit_line,'-k')
%     
%     if ismember(icondition,[1,3])
%         ylabel('Hit Rate (%)','FontSize',18)
%     end
%     if ismember(icondition,[3,4])
%         xlabel('D-prime','FontSize',18)
%     end
%     title(leg_labels(icondition),'FontSize',18)
%     xlim([-2,5.5])
%     ylim([0,1])
%     text(1,0.2,"R^2 = " + num2str(R_squared))
%     text(1,0.1,"P = " + num2str(P(2)))
% 
% end

