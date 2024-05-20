%% srm_nirs_eeg_compare_beta_with_behavior

% load behavior data
load('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_Behavior_Results.mat')
% ORDER: itd50, itd500, ildnat, ild10

% load beta values
load('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\collapsed attend and masker PFC time constant SPEECH NO BREATH STG results.mat')
stg_betas_to_plot = all_betas;
stg_means_to_plot = all_means;

load('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\collapsed attend and masker PFC time constant SPEECH NO BREATH PFC results.mat')
pfc_betas_to_plot = all_betas;
pfc_means_to_plot = all_means;

% ORDER: itd50, itd500, ild10, ildnat

% Define current subjects
subject_ID = char('NDARVX753BR6','NDARZD647HJ1','NDARBL382XK5','NDARGF569BF3','NDARBA306US5','NDARFD284ZP3','NDARAS648DT4','NDARLM531OY3','NDARXL287BE1','NDARRF358KO3','NDARGT639XS6','NDARDC882NK4','NDARWB491KR3','NDARNL224RR9','NDARTT639AB1','NDARAZC45TW3','NDARNS784LM2','NDARLB144ZM4','NDARTP382XC8','NDARLJ581GD7','NDARGS283RM9'); %

conditions = string({'itd50','itd500','ildnat','ild10'});
colors = {'r','g','b','m'};
%for icondition = 1:length(conditions)
    for isubject = 1:length(subject_ID)
        % get d-prime for this subject and this condition
        this_subject_d_prime = mean(d_primes_speech_masker(:,isubject),1);

        % get pfc beta for this subject and this condition
        this_subject_pfc_beta = mean(pfc_betas_to_plot(isubject,:),2);
        % get stg beta for this subject and this condition
        this_subject_stg_beta = mean(stg_betas_to_plot(isubject,:),2);
        % get pfc AUC for this subject and this condition
        this_subject_pfc_mean = mean(pfc_means_to_plot(isubject,:),2);
        % get pfc AUC for this subject and this condition
        this_subject_stg_mean = mean(stg_means_to_plot(isubject,:),2);

        figure(1) % PFC Beta
        scatter(this_subject_d_prime,this_subject_pfc_beta,'filled','MarkerFaceColor','m')
        hold on
        figure(2) % stg Beta
        scatter(this_subject_d_prime,this_subject_stg_beta,'filled','MarkerFaceColor','m')
        hold on
        figure(3) % pfc AUC
        scatter(this_subject_d_prime,this_subject_pfc_mean,'filled','MarkerFaceColor','m')
        hold on
        figure(4) % stg AUC
        scatter(this_subject_d_prime,this_subject_stg_mean,'filled','MarkerFaceColor','m')
        hold on
    end
%end
figure(1)
% best fit lines
%for icondition = 1:length(conditions)
    x = mean(d_primes_speech_masker,1);%(icondition,:);

    y = mean(pfc_betas_to_plot,2);%(:,icondition);

    coefficients = polyfit(x, y, 1);
    % Get the estimated yFit value for each of those 1000 new x locations.
    yFit = polyval(coefficients , x)';
    p = plot(x,yFit,'m');
%end
% Calculate R^2 value
y_mean = mean(y);
SS_tot = sum((y - y_mean).^2);
SS_res = sum((y - yFit).^2);
R_squared = 1 - SS_res / SS_tot;
text(2,0.1,['R2 = ',string(R_squared)])
% plot labels
xlabel("Behavioral Sensitivity (d')",'FontSize',18)
ylabel('PFC Beta','FontSize',18)
title('PFC Beta vs. D-prime','FontSize',18)

figure(2)
% best fit lines
%for icondition = 1:length(conditions)
    x = mean(d_primes_speech_masker,1);%(icondition,:);

    y = mean(stg_betas_to_plot,2);%(:,icondition);

    coefficients = polyfit(x, y, 1);

    % Get the estimated yFit value for each of those 1000 new x locations.
    yFit = polyval(coefficients , x)';
    plot(x,yFit,'m')
%end

% Calculate R^2 value
y_mean = mean(y);
SS_tot = sum((y - y_mean).^2);
SS_res = sum((y - yFit).^2);
R_squared = 1 - SS_res / SS_tot;
text(2,0.1,['R2 = ',string(R_squared)])
%plot labels
xlabel("Behavioral Sensitivity (d')",'FontSize',18)
ylabel('STG Beta','FontSize',18)
title('STG Beta vs. D-prime','FontSize',18)

figure(3)
% best fit lines
%for icondition = 1:length(conditions)
    x = mean(d_primes_speech_masker,1);%(icondition,:);

    y = mean(pfc_means_to_plot,2);%(:,icondition);

    coefficients = polyfit(x, y, 1);

    % Get the estimated yFit value for each of those 1000 new x locations.
    yFit = polyval(coefficients , x)';
    plot(x,yFit,'m')
%end
% Calculate R^2 value
y_mean = mean(y);
SS_tot = sum((y - y_mean).^2);
SS_res = sum((y - yFit).^2);
R_squared = 1 - SS_res / SS_tot;
text(2,0.1,['R2 = ',string(R_squared)])
%plot labels
xlabel("Behavioral Sensitivity (d')",'FontSize',18)
ylabel('PFC Mean','FontSize',18)
title('PFC Mean vs. D-prime','FontSize',18)


figure(4)
% best fit lines
%for icondition = 1:length(conditions)
    x = mean(d_primes_speech_masker,1); %(icondition,:);

    y = mean(stg_means_to_plot,2); %(:,icondition);

    coefficients = polyfit(x, y, 1);
    % Get the estimated yFit value for each of those 1000 new x locations.
    yFit = polyval(coefficients , x)';
    plot(x,yFit,'m')
%end
% Calculate R^2 value
y_mean = mean(y);
SS_tot = sum((y - y_mean).^2);
SS_res = sum((y - yFit).^2);
R_squared = 1 - SS_res / SS_tot;
text(2,0.1,['R2 = ',string(R_squared)])
%plot labels
xlabel("Behavioral Sensitivity (d')",'FontSize',18)
ylabel('STG Mean','FontSize',18)
title('STG Mean vs. D-prime','FontSize',18)
