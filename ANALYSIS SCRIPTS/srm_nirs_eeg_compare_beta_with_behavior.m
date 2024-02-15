%% srm_nirs_eeg_compare_beta_with_behavior

% load behavior data
load('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\SRM-NIRS-EEG-1_Behavior_Results.mat')
% ORDER: itd50, itd500, ildnat, ild10

% load beta values
load('C:\Users\benri\Documents\GitHub\SRM-NIRS-EEG\RESULTS DATA\collapsed attend and masker PFC time constant SPEECH results.mat')
% ORDER: control, itd50, itd500, ild10, ildnat
pfc_betas_to_plot(:,1) = [];
stg_betas_to_plot(:,1) = [];
all_AUC_stg(:,1) = [];
all_AUC_pfc(:,1) = [];

% Define current subjects
subject_ID = char('NDARVX753BR6','NDARZD647HJ1','NDARBL382XK5','NDARGF569BF3','NDARBA306US5','NDARFD284ZP3','NDARAS648DT4','NDARLM531OY3','NDARXL287BE1','NDARRF358KO3','NDARGT639XS6','NDARDC882NK4','NDARWB491KR3','NDARNL224RR9','NDARTT639AB1','NDARAZC45TW3','NDARNS784LM2','NDARLB144ZM4','NDARTP382XC8','NDARLJ581GD7','NDARGS283RM9'); %

conditions = string({'itd50','itd500','ildnat','ild10'});
colors = {'r','g','b','m'};
for icondition = 1:length(conditions)
    for isubject = 1:length(subject_ID)
        % get d-prime for this subject and this condition
        this_subject_d_prime = d_primes_speech_masker(icondition,isubject);

        % get pfc beta for this subject and this condition
        this_subject_pfc_beta = pfc_betas_to_plot(isubject,icondition);
        % get stg beta for this subject and this condition
        this_subject_stg_beta = stg_betas_to_plot(isubject,icondition);
        % get pfc AUC for this subject and this condition
        this_subject_pfc_AUC = all_AUC_pfc(isubject,icondition);
        % get pfc AUC for this subject and this condition
        this_subject_stg_AUC = all_AUC_stg(isubject,icondition);

        figure(1) % PFC Beta
        scatter(this_subject_d_prime,this_subject_pfc_beta,'filled','MarkerFaceColor',string(colors(icondition)))
        hold on
        figure(2) % stg Beta
        scatter(this_subject_d_prime,this_subject_stg_beta,'filled','MarkerFaceColor',string(colors(icondition)))
        hold on
        figure(3) % pfc AUC
        scatter(this_subject_d_prime,this_subject_pfc_AUC,'filled','MarkerFaceColor',string(colors(icondition)))
        hold on
        figure(4) % stg AUC
        scatter(this_subject_d_prime,this_subject_stg_AUC,'filled','MarkerFaceColor',string(colors(icondition)))
        hold on
    end
end
figure(1)
% best fit lines
for icondition = 1:length(conditions)
    x = d_primes_speech_masker(icondition,:);

    y = pfc_betas_to_plot(:,icondition);

    coefficients = polyfit(x, y, 1);
    % Create a new x axis with exactly 1000 points (or whatever you want).
    xFit = linspace(min(x), max(x), 1000);
    % Get the estimated yFit value for each of those 1000 new x locations.
    yFit = polyval(coefficients , xFit);
    p(icondition,:) = plot(xFit,yFit,string(colors(icondition)));
end
% plot labels
xlabel("Behavioral Sensitivity (d')",'FontSize',18)
ylabel('PFC Beta','FontSize',18)
title('PFC Beta vs. D-prime','FontSize',18)
legend([p(1),p(2),p(3),p(4)],conditions)

figure(2)
% best fit lines
for icondition = 1:length(conditions)
    x = d_primes_speech_masker(icondition,:);

    y = stg_betas_to_plot(:,icondition);

    coefficients = polyfit(x, y, 1);
    % Create a new x axis with exactly 1000 points (or whatever you want).
    xFit = linspace(min(x), max(x), 1000);
    % Get the estimated yFit value for each of those 1000 new x locations.
    yFit = polyval(coefficients , xFit);
    plot(xFit,yFit,string(colors(icondition)))
end
%plot labels
xlabel("Behavioral Sensitivity (d')",'FontSize',18)
ylabel('STG Beta','FontSize',18)
title('STG Beta vs. D-prime','FontSize',18)
legend([p(1),p(2),p(3),p(4)],conditions)

figure(3)
% best fit lines
for icondition = 1:length(conditions)
    x = d_primes_speech_masker(icondition,:);

    y = all_AUC_pfc(:,icondition);

    coefficients = polyfit(x, y, 1);
    % Create a new x axis with exactly 1000 points (or whatever you want).
    xFit = linspace(min(x), max(x), 1000);
    % Get the estimated yFit value for each of those 1000 new x locations.
    yFit = polyval(coefficients , xFit);
    plot(xFit,yFit,string(colors(icondition)))
end
%plot labels
xlabel("Behavioral Sensitivity (d')",'FontSize',18)
ylabel('PFC AUC','FontSize',18)
title('PFC AUC vs. D-prime','FontSize',18)
legend([p(1),p(2),p(3),p(4)],conditions)


figure(4)
% best fit lines
for icondition = 1:length(conditions)
    x = d_primes_speech_masker(icondition,:);

    y = all_AUC_stg(:,icondition);

    coefficients = polyfit(x, y, 1);
    % Create a new x axis with exactly 1000 points (or whatever you want).
    xFit = linspace(min(x), max(x), 1000);
    % Get the estimated yFit value for each of those 1000 new x locations.
    yFit = polyval(coefficients , xFit);
    plot(xFit,yFit,string(colors(icondition)))
end
%plot labels
xlabel("Behavioral Sensitivity (d')",'FontSize',18)
ylabel('STG AUC','FontSize',18)
title('STG AUC vs. D-prime','FontSize',18)
legend([p(1),p(2),p(3),p(4)],conditions)
