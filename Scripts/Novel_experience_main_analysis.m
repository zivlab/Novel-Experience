clear;clc;close all
%% Load processed data and settings:

% Load processed data:
data_path = 'D:\Gal PhD\Projects\Experience_Project\scripts final github upload\processed data\'; %change to downloaded data directory
cd(data_path)

load([data_path,'example_cells_data.mat'])
load([data_path,'PV_correlation_mats_data.mat'])
load([data_path,'Comparing_corr_across_phases_data.mat'])
load([data_path,'Parameters_across_days_data.mat'])
load([data_path,'Comparing_corr_across_phases_bygroup_data.mat'])
load([data_path,'Comparing_corr_across_phases_byphase_data.mat'])
load([data_path,'Comparing_across_phases_ensemble_overlap_data.mat'])
load([data_path,'Laps_control_original_vs_subsampaled_data.mat'])
load([data_path,'Comparing_corr_across_phases_lap_control_data.mat'])

% Setting parameters:
phase1 = 1:6;
phase2 = 6:11;
phase3 = 11:16;

phases_labels_name = {'Baseline' 'Novel env. 15 min' 'Novel env. 24 hrs'};
phases_labels = {'Phase 1' 'Phase 2' 'Phase 3'};

corr_var_label = {'Tuning Curve Correlation','Ensemble Rate Correlation'};

%% Figure 1C - Example place cells:

ylims = [7402.9  9887.8 ; 6118.5  8000; 10231 11858;1597  3597; 2229.6  4229.6; 709  2067.2;...
1327.4    3327.4; 5466.7    7066.7; 4363.4    5736.6; 8819    10685; 5744.4    7211.1; 5160.7    6354.8];
for cell = 1:4
    f = figure;
    f.Position =  [795 144 216 623];
    for day = 1:3
        subplot(3,1,day)
        location = example_cells_data.cells(cell).location{day};
        active_position = example_cells_data.cells(cell).activity{day};
        plot(location,1:length(location),'Color','[0.8 0.8 0.8]')
        hold on
        plot(location(active_position),active_position,'.','MarkerSize', 10,'Color','[0.0157    0.2431    0.5804]')
        xlim([0 26])
        ylim([ylims(day+3*(cell-1),1) ylims(day+3*(cell-1),2)])
        yticklabels([])
        xticklabels([])
        ylabel(phases_labels_name(day))
        if day==1
            title(['cell ' num2str(cell)])
        end
        box off
        axis square
    end
end

%% Figure 1D - Population vector correlations in spatial representation across recording days:

% Upper panel - Population vector correlations between positions on the linear track across 16 recording days: 

bins = 20;
sessions_line_diag = bins+0.5:bins:(bins+0.5)*15;
PV_mean_all = median(PV_corr_mice,3,'omitnan');
f = figure;
f.Position = [2297 78 933 810];
imagesc(PV_mean_all,[0 0.73]) 
colormap(jet)
colorbar
xline(sessions_line_diag,'LineWidth',1.5,'Color','k','Alpha',1);
yline(sessions_line_diag,'LineWidth',1.5,'Color','k','Alpha',1);
xticks(10:20:350)
xticklabels(1:2:31)
xlabel('Days')
yticks(10:20:350)
yticklabels(1:2:31)
ylabel('Days')
ytickangle(90)
axis square

% Lower-left panel - Population vector correlations across recording days averaged across matching positions: 

PV_mean_phases_mean_all = median(PV_corr_mean_mice,3,'omitnan');
figure
imagesc(PV_mean_phases_mean_all,[0.35 0.71]) 
colormap(jet)
colorbar
xticks(1:16)
xlabel('Days')
yticks(1:16)
ylabel('Days')
ytickangle(90)
axis square

% Lower-right panel - Population vector correlations across positions and running directions averaged across all pairs of sessions: 

PV_corr_both_all = median(PV_corr_both_mean_mice,3,'omitnan');
figure
imagesc(PV_corr_both_all,[0 0.55])
colormap(jet)
colorbar
xline(20.5,'LineWidth',3,'Color','w','Alpha',1);
yline(20.5,'LineWidth',3,'Color','w','Alpha',1);
xticks(10:20:50)
xticklabels({'Rightward','Leftward'})
yticks(10:20:50)
yticklabels({'Rightward','Leftward'})
ytickangle(90)
axis square

%%  Figure 2A-B, Right - Comparing drift in spatial tuning and activity rates between baseline (phase 1) and novel enviornment visits (phase 2&3):
clc
for c = 2:3
    figure
    drift_phase1 = squeeze(drift_bs_vs_enr(1,:,c,:));
    drift_phase1_within = [mean(within_all(:,phase1,c),2,'omitnan'),drift_phase1'];
    std_drift_phase1 = std(drift_phase1_within,0,1,'omitnan')/sqrt(8);
    errorbar(mean(drift_phase1_within,'omitnan'),std_drift_phase1,'LineWidth',1.5,'Color',[0.65098    0.65098    0.65098]);
    hold on
    drift_phase2_3 = squeeze(drift_bs_vs_enr(2,:,c,:));
    drift_phase2_3_within_mean = mean([mean(within_all(:,phase2,c),2,'omitnan'),mean(within_all(:,phase3,c),2,'omitnan')],2);
    drift_phase2_3_within = [drift_phase2_3_within_mean,drift_phase2_3'];
    std_drift_phase2_3 = std(drift_phase2_3_within,0,1,'omitnan')/sqrt(8);
    errorbar(mean(drift_phase2_3_within,'omitnan'),std_drift_phase2_3,'LineWidth',1.5,'Color',[0.72941      0.24706      0.24706]);
    legend([{'Baseline'},{'Novel env'}])
    legend boxoff
    if c==3
        ylim([0.45 1])
    else
        ylim([0.35 0.85])
    end
    ylabel(corr_var_label{c-1})
    xlim([0 7])
    xticks(1:6)
    xticklabels(0:2:10)
    xlabel('Elapsed time (days)')
    box off
    axis square

   % check sig bw manipulation:
   bs_m = drift_phase1_within';
   enr_m = drift_phase2_3_within';
   [p,tbl,~] = friedman([bs_m(:) enr_m(:)],6);
   p_val_bw_manipulation = round(p,3);
   chi_bw_manipulation = round(tbl{2,5},3);

   disp(['Comparing drift in ' corr_var_label{c-1} ' of baseline vs visits to novel environments'])
   disp(['Using the nonparametric Friedman test: p-value = ' num2str(p_val_bw_manipulation) ' Chi = ' num2str(chi_bw_manipulation)])
end

%% Figure 2A-B, Left - Comparing spatial tuning and activity rates across all pairs of sessions regardless of elapsed time between baseline (phase 1) and novel enviornment visits (phase 2&3):

mean_across = mean(mean_all_sessions,2,'omitnan');

clc
for c = 2:3
    mean_phase1 = squeeze(mean_across(1,:,c,:));
    mean_phase2_3 = mean([squeeze(mean_across(2,:,c,:)) ,squeeze(mean_across(3,:,c,:))],2,'omitnan');
    bs_phases = [mean_phase1,mean_phase2_3];
    std_across_mice = std(bs_phases,'omitnan')./sqrt(8);
    figure
    hold on
    plot(bs_phases','color',[0.6 0.6 0.6])
    errorbar(mean(bs_phases,'omitnan'),std_across_mice,'color','k','LineWidth',1.5)
    if c==3
        ylim([0.45 1])
    else
        ylim([0.35 0.85])
    end
    xlim([0.5 2.5])
    ylabel(corr_var_label{c-1})
    set(gca,'xtick',1:2,'xticklabel',[{'Baseline'},{'Novel env'}])
    box off
    axis square

    % check sig wilcoxson:
    [p_val,~,~] = ranksum(bs_phases(:,1),bs_phases(:,2));
    p_val_wilcxon = round(p_val,3);

    disp(['Comparing ' corr_var_label{c-1} ' across all pairs of sessions regardless of elapsed time of baseline vs visits to novel environments'])
    disp(['Using the nonparametric Wilcoxon signed-rank test: p-value = ' num2str(p_val_wilcxon)])
end

%% Figure S1A, Right - Comparing drift in population vector activity between baseline (phase 1) and novel enviornment visits (phase 2&3):

c = 1;
figure
drift_phase1 = squeeze(drift_bs_vs_enr(1,:,c,:));
drift_phase1_within = [mean(within_all(:,phase1,c),2,'omitnan'),drift_phase1'];
std_drift_phase1 = std(drift_phase1_within,0,1,'omitnan')/sqrt(8);
errorbar(mean(drift_phase1_within,'omitnan'),std_drift_phase1,'LineWidth',1.5,'Color',[0.65098    0.65098    0.65098]);
hold on
drift_phase2_3 = squeeze(drift_bs_vs_enr(2,:,c,:));
drift_phase2_3_within_mean = mean([mean(within_all(:,phase2,c),2,'omitnan'),mean(within_all(:,phase3,c),2,'omitnan')],2);
drift_phase2_3_within = [drift_phase2_3_within_mean,drift_phase2_3'];
std_drift_phase2_3 = std(drift_phase2_3_within,0,1,'omitnan')/sqrt(8);
errorbar(mean(drift_phase2_3_within,'omitnan'),std_drift_phase2_3,'LineWidth',1.5,'Color',[0.72941      0.24706      0.24706]);
legend([{'Baseline'},{'Novel env'}])
legend boxoff
ylim([0.35 0.85])
ylabel('Population vector correlation')
xlim([0 7])
xticks(1:6)
xticklabels(0:2:10)
xlabel('Elapsed time (days)')
box off
axis square

% check sig bw manipulation:
bs_m = drift_phase1_within';
enr_m = drift_phase2_3_within';
[p,tbl,~] = friedman([bs_m(:) enr_m(:)],6);
p_val_bw_manipulation = round(p,3);
chi_bw_manipulation = round(tbl{2,5},3);

clc
disp('Comparing drift in Population vector correlation of baseline vs visits to novel environments')
disp(['Using the nonparametric Friedman test: p-value = ' num2str(p_val_bw_manipulation) ' Chi = ' num2str(chi_bw_manipulation)])

%% Figure S1A, Left - Comparing  population vector acticity across all pairs of sessions regardless of elapsed time between baseline (phase 1) and novel enviornment visits (phase 2&3):

mean_across = mean(mean_all_sessions,2,'omitnan');
c = 1;
mean_phase1 = squeeze(mean_across(1,:,c,:));
mean_phase2_3 = mean([squeeze(mean_across(2,:,c,:)) ,squeeze(mean_across(3,:,c,:))],2,'omitnan');
bs_phases = [mean_phase1,mean_phase2_3];
std_across_mice = std(bs_phases,'omitnan')./sqrt(8);
figure
hold on
plot(bs_phases','color',[0.6 0.6 0.6])
errorbar(mean(bs_phases,'omitnan'),std_across_mice,'color','k','LineWidth',1.5)
ylim([0.35 0.85])
xlim([0.5 2.5])
ylabel('Population vector correlation')
set(gca,'xtick',1:2,'xticklabel',[{'Baseline'},{'Novel env'}])
box off
axis square

% check sig wilcoxson:
[p_val,~,~] = ranksum(bs_phases(:,1),bs_phases(:,2));
p_val_wilcxon = round(p_val,3);

clc
disp('Comparing Population vector correlation across all pairs of sessions regardless of elapsed time of baseline vs visits to novel environments')
disp(['Using the nonparametric Wilcoxon signed-rank test: p-value = ' num2str(p_val_wilcxon)])
 
%% Figure S1B-E - Gross statistical attributes of the spatial code across recording days:

% Panel B - number of place cells
% Panel C - mean event rate
% Panel D - spatial information
% Panel E - within-day stability (even vs odd laps correlation)

parameters = {mean_activity_mice,sp_spike_mice,num_sig_cells_mice,pv_within_lap};
ylabels = {'Mean activity (mean events)','Mean spatial information (bit/event)','Place cells','Within day correlation'};
ylims = { [0 0.005];[0 3.2] ; [0 400] ; [0 1]};

clc
for c = 1:length(parameters)
    var_curr = parameters{c};
    figure
    errorbar(mean(var_curr,'omitnan'),std(var_curr,'omitnan')/sqrt(8),'Color','k','LineWidth',2)
    xticks(1:16)
    xlim([0 17])
    xline([6 11],'LineWidth',1.5,'Alpha',0.2)
    ylabel(ylabels{c})
    ylim(ylims{c})
    box off
    axis square

    % check sig wilcoxson:
    phase1_mean = mean(var_curr(:,phase1),2,'omitnan');
    phase2_mean = mean(var_curr(:,phase2),2,'omitnan');
    phase3_mean = mean(var_curr(:,phase3),2,'omitnan');
    all_phase_mean = [phase1_mean phase2_mean phase3_mean];
    [p_val_phase1_2,~,~] = ranksum(all_phase_mean(:,1),all_phase_mean(:,2));
    [p_val_phase1_3,~,~] = ranksum(all_phase_mean(:,1),all_phase_mean(:,3));
    [p_val_phase2_3,~,~] = ranksum(all_phase_mean(:,2),all_phase_mean(:,3));
   
    disp([ylabels{c} ' across phases'])
    disp(['Using the nonparametric Wilcoxon signed-rank test: p-value b/w phase 1 and 2 = ' num2str(round(p_val_phase1_2,3))])
    disp(['Using the nonparametric Wilcoxon signed-rank test: p-value b/w phase 1 and 3 = ' num2str(round(p_val_phase1_3,3))])
    disp(['Using the nonparametric Wilcoxon signed-rank test: p-value b/w phase 2 and 3 = ' num2str(round(p_val_phase2_3,3))])
end

%% Figure S2A-D, Right panels - Comparing drift between baseline (phase 1) and novel enviornment visits separately for each phase as a function of elapsed time:
% Figure S2A, Right panel - Comparing phase 1 with 15 minutes after novel env tuning curve correlations 
% Figure S2B, Right panel - Comparing phase 1 with 15 minutes after novel env ensemble rate correlations
% Figure S2C, Right panel - Comparing phase 1 with 24 hours after novel env tuning curve correlations 
% Figure S2D, Right panel - Comparing phase 1 with 24 hours after novel env ensemble rate correlations 

phases = {phase1,phase2,phase3};

clc
for p = 1:2
    for c = 2:3
        figure
        drift_phase1 = squeeze(drift_bs_vs_enr_bygroup(1,:,c,:));
        drift_phase1_within = [mean(within_all_bygroup(:,phases{1},c),2,'omitnan'),drift_phase1'];
        std_drift_phase1 = std(drift_phase1_within,0,1,'omitnan')/sqrt(8);
        errorbar(mean(drift_phase1_within,'omitnan'),std_drift_phase1,'LineWidth',1.5,'Color',[0.65098    0.65098    0.65098]);
        hold on
        drift_phase2_3_bygroup = squeeze(drift_bs_vs_enr_bygroup(p+1,:,c,:));
        drift_phase2_3_within_mean_bygroup = mean(within_all_bygroup(:,phases{p+1},c),2,'omitnan');
        drift_phase2_3_within_bygroup = [drift_phase2_3_within_mean_bygroup,drift_phase2_3_bygroup'];
        std_drift_phase2_3_bygroup = std(drift_phase2_3_within_bygroup,0,1,'omitnan')/sqrt(8);
        errorbar(mean(drift_phase2_3_within_bygroup,'omitnan'),std_drift_phase2_3_bygroup,'LineWidth',1.5,'Color',[0.72941      0.24706      0.24706]);
        legend([{'Baseline'},phases_labels_name{p+1}])
        legend boxoff
        if c==3
            ylim([0.45 1])
        else
            ylim([0.35 0.85])
        end
        ylabel(corr_var_label{c-1})
        xlim([0 7])
        xticks(1:6)
        xticklabels(0:2:10)
        xlabel('Elapsed time (days)')
        box off
        axis square

        % check sig:
        bs_m = drift_phase1_within';
        enr_m = drift_phase2_3_within_bygroup';
        [pval,tbl,~] = friedman([bs_m(:) enr_m(:)],6);
        p_val_bw_manipulation = round(pval,3);
        chi_bw_manipulation = round(tbl{2,5},3);

        disp(['Comparing drift in ' corr_var_label{c-1} ' of baseline vs exposure to ' phases_labels_name{p+1}])
        disp(['Using the nonparametric Friedman test: p-value = ' num2str(p_val_bw_manipulation) ' Chi = ' num2str(chi_bw_manipulation)])
    end
end

%% Figure S2A-D, Left panels -  Comparing correlations across all pairs of sessions regardless of elapsed time between baseline (phase 1) and novel enviornment visits separately for each phase:
% Figure S2A, Left panel - Comparing phase 1 with 15 minutes after novel env tuning curve correlations 
% Figure S2B, Left panel - Comparing phase 1 with 15 minutes after novel env ensemble rate correlations 
% Figure S2C, Left panel - Comparing phase 1 with 24 hours after novel env tuning curve correlations 
% Figure S2D, Left panel - Comparing phase 1 with 24 hours after novel env ensemble rate correlations 

mean_across = mean(mean_all_sessions_bygroup,2,'omitnan');

clc
for p = 1:2
    for c = 2:3
        mean_phase1 = squeeze(mean_across(1,:,c,:));
        mean_phase_phase2_3 = squeeze(mean_across(p+1,:,c,:));
        bs_phases = [mean_phase1, mean_phase_phase2_3];
        std_across_mice = std(bs_phases,'omitnan')./sqrt(8);
        figure
        hold on
        plot(bs_phases','color',[0.6 0.6 0.6])
        errorbar(mean(bs_phases,'omitnan'),std_across_mice,'color','k','LineWidth',1.5)
        if c==3
            ylim([0.45 1])
        else
            ylim([0.35 0.85])
        end
        xlim([0.5 2.5])
        ylabel(corr_var_label{c-1})
        set(gca,'xtick',1:2,'xticklabel',[{'Baseline'},phases_labels_name{p+1}])
        box off
        axis square

        % check sig wilcoxson:
        [p_val,~,~] = ranksum(bs_phases(:,1),bs_phases(:,2));
        p_val_wilcxon = round(p_val,3);

        disp(['Comparing ' corr_var_label{c-1} ' across all pairs of sessions regardless of elapsed time of baseline vs exposure to ' phases_labels_name{p+1}])
        disp(['Using the nonparametric Wilcoxon signed-rank test: p-value = ' num2str(p_val_wilcxon)])
    end
end

%% Figure S2E - Comparing novel exposure scheduals across all pairs of sessions regardless of elapsed time:

% Right panel - Comparing novel exposure after 15 minutes with 24 hours tuning curve correlations
% Left panel - Comparing novel exposure after 15 minutes with 24 hours ensemble rate correlations

mean_across = mean(mean_all_sessions_bygroup,2,'omitnan');

clc
for c = 2:3
    mean_15min = squeeze(mean_across(2,:,c,:));
    mean_24hrs = squeeze(mean_across(3,:,c,:));
    bs_phases = [mean_15min,mean_24hrs];
    std_across_mice = std(bs_phases,'omitnan')./sqrt(8);
    figure
    hold on
    plot(bs_phases','color',[0.6 0.6 0.6])
    errorbar(mean(bs_phases,'omitnan'),std_across_mice,'color','k','LineWidth',1.5)
    if c==3
        ylim([0.45 1])
    else
        ylim([0.35 0.85])
    end
    xlim([0.5 2.5])
    ylabel(corr_var_label{c-1})
    set(gca,'xtick',1:2,'xticklabel',[phases_labels_name(2),phases_labels_name(3)])
    box off
    axis square

    % check sig wilcoxson:
    [p_val,~,~] = ranksum(bs_phases(:,1),bs_phases(:,2));
    p_val_wilcxon_control = round(p_val,3);

    disp(['Comparing ' corr_var_label{c-1} ' across all pairs of sessions regardless of elapsed time of exposure to novel env 15 mins vs 24 hours after'])
    disp(['Using the nonparametric Wilcoxon signed-rank test: p-value = ' num2str(p_val_wilcxon_control)])
end

%% Figure S2F - Comparing phases 2 & 3 of experiment across all pairs of sessions regardless of elapsed time:

% Right panel - Comparing phase 2 with 3 tuning curve correlations 
% Left panel - Comparing phase 2 with 3 ensemble rate correlations

mean_across = mean(mean_all_sessions_byphase,2,'omitnan');

clc
for c = 1:2
    mean_phase2 = squeeze(mean_across(2,:,c,:));
    mean_phase3 = squeeze(mean_across(3,:,c,:));
    bs_phases = [mean_phase2, mean_phase3];
    std_across_mice = std(bs_phases,'omitnan')./sqrt(8);
    figure
    hold on
    plot(bs_phases','color',[0.6 0.6 0.6])
    errorbar(mean(bs_phases,'omitnan'),std_across_mice,'color','k','LineWidth',1.5)
    if c==2
        ylim([0.45 1])
    else
        ylim([0.35 0.85])
    end
    xlim([0.5 2.5])
    ylabel(corr_var_label{c})
    set(gca,'xtick',1:2,'xticklabel',[phases_labels(2),phases_labels(3)])
    box off
    axis square

    % check sig wilcoxson:
    [p_val,~,~] = ranksum(bs_phases(:,1),bs_phases(:,2));
    p_val_wilcxon_control = round(p_val,3);

    disp(['Comparing ' corr_var_label{c} ' across all pairs of sessions regardless of elapsed time of phase 2 vs 3'])
    disp(['Using the nonparametric Wilcoxon signed-rank test: p-value = ' num2str(p_val_wilcxon_control)])
end

%% Figure S2G - Comparing overlap percentage of ensmble activity:

% Right panel - Comparing phase 1 with 2&3 overlap percentage as a function of elapsed time: 

drift_phase1_R = squeeze(drift_bs_vs_enr_overlap(1,:,1,:));
drift_phase1_L = squeeze(drift_bs_vs_enr_overlap(1,:,2,:));
figure
std_mat_bs = [drift_phase1_R drift_phase1_L]*100;
std_phase = std(std_mat_bs,0,2,'omitnan')/sqrt(8);
errorbar(mean(std_mat_bs,2,'omitnan'),std_phase,'LineWidth',1.5,'Color',[0.65098    0.65098    0.65098]);
hold on
drift_phase2_3_R = squeeze(drift_bs_vs_enr_overlap(2,:,1,:));
drift_phase2_3_L = squeeze(drift_bs_vs_enr_overlap(2,:,2,:));
std_mat_p = [drift_phase2_3_R drift_phase2_3_L]*100;
std_phase = std(std_mat_p,0,2,'omitnan')/sqrt(8);
errorbar(mean(std_mat_p,2,'omitnan'),std_phase,'LineWidth',1.5,'Color',[0.72941      0.24706      0.24706]);
legend([{'Baseline'},{'Novel env'}])
legend boxoff
ylim([10 60])
ylabel('Percentage overlap')
xlim([0 6])
xticks(1:5)
xticklabels(2:2:10)
xlabel('Elapsed time (days)')
box off
axis square

% check sig bw manipulation:

bs_m = std_mat_bs';
enr_m = std_mat_p';
[p,tbl,stats] = friedman([bs_m(:) enr_m(:)],10);
p_val_bw_manipulation = round(p,3);
chi_bw_manipulation = round(tbl{2,5},3);

clc
disp('Comparing ensemble overlap percentage of baseline vs visits to novel environments')
disp(['Using the nonparametric Friedman test: p-value = ' num2str(p_val_bw_manipulation) ' Chi = ' num2str(chi_bw_manipulation)])

% Left panel - Comparing phase 1 with 2&3 overlap percentages across all pairs of sessions regardless of elapsed time: 

mean_across = mean(mean_all_sessions_overlap,2,'omitnan');
mean_across_sides = mean(mean_across,3,'omitnan');
mean_phase1 = squeeze(mean_across_sides(1,:,1,:));
mean_phase2_3 = mean([squeeze(mean_across_sides(2,:,1,:)) ,squeeze(mean_across(3,:,1,:))],2,'omitnan');
bs_phases = [mean_phase1, mean_phase2_3]*100;
std_across_mice = std(bs_phases,'omitnan')./sqrt(8);
figure
hold on
plot(bs_phases','color',[0.6 0.6 0.6])
errorbar(mean(bs_phases,'omitnan'),std_across_mice,'color','k','LineWidth',1.5)
ylim([10 60])
xlim([0.5 2.5])
ylabel('Percentage overlap')
set(gca,'xtick',1:2,'xticklabel',[{'Baseline'},{'Novel env'}])
box off
axis square

% check sig wilcoxson:
[p_val,~,~] = ranksum(bs_phases(:,1),bs_phases(:,2));
p_val_wilcxon_control = round(p_val,3);

disp('Comparing ensemble overlap percentage across all pairs of sessions regardless of elapsed time of baseline vs visits to novel environments')
disp(['Using the nonparametric Wilcoxon signed-rank test: p-value = ' num2str(p_val_wilcxon_control)])

%% Figure S3A - The mean number of laps transversed per mouse in each phase in the original and subsampled data:

% Upper panel - Original data 
% Lower panel - Subsampled data

num_laps = {num_laps_all  num_laps_chosen};
for comp = 1:2
    figure
    num_laps_all_mice = num_laps{comp};
    baseline_mean_num_laps = mean(num_laps_all_mice(:,1:6),2);
    p1_mean = mean(num_laps_all_mice(:,7:11),2);
    p2_mean = mean(num_laps_all_mice(:,12:16),2);
    all_phase_mean = [baseline_mean_num_laps p1_mean p2_mean]; 
    plot(all_phase_mean','color',[0.6 0.6 0.6])
    hold on
    errorbar(mean(all_phase_mean),std(all_phase_mean)/sqrt(8),'color','k','LineWidth',1.5)
    xlim([0.2 3.8])
    xticks(1:3)
    y_min = min(all_phase_mean,[],'all');
    y_max = max(all_phase_mean,[],'all');
    ylim([0  300])
    xticklabels(phases_labels)
    ylabel('Laps')
    box off
    axis square

end

% change percent:
phase_both = 7:16;
baseline_mean_num_laps = mean(num_laps_all(:,phase1),2);
p2_3_mean = mean(num_laps_all(:,phase_both),2);
per_change_laps = p2_3_mean./baseline_mean_num_laps;
std_laps = std(per_change_laps)/sqrt(8);
mean_per_change_laps = mean(per_change_laps);

clc
disp(['Percent change in averaged lap number between phase 1 and phase 2&3: ' num2str(round(mean_per_change_laps*10,0)) ' ± ' num2str(round(std_laps,2)) '%'])

%% Figure S3B-C, Right - Comparing drift in spatial tuning and activity rates in subsampled data (same number of laps across days):
 
% Upper panel - Comparing phase 1 with 2&3 tuning curve correlations as a function of elapsed time
% Lower panel - Comparing phase 1 with 2&3 ensemble rate correlations as a function of elapsed time 

clc
for c = 1:2
    figure
    drift_phase1 = squeeze(drift_bs_vs_enr_lap_control(1,:,c,:));
    drift_phase1_within = [mean(within_all_lap_control(:,phase1,c),2,'omitnan'),drift_phase1'];
    std_drift_phase1 = std(drift_phase1_within,0,1,'omitnan')/sqrt(8);
    errorbar(mean(drift_phase1_within,'omitnan'),std_drift_phase1,'LineWidth',1.5,'Color',[0.65098    0.65098    0.65098]);
    hold on
    drift_phase2_3 = squeeze(drift_bs_vs_enr_lap_control(2,:,c,:));
    drift_phase2_3_within_mean = mean([mean(within_all_lap_control(:,phase2,c),2,'omitnan'),mean(within_all_lap_control(:,phase3,c),2,'omitnan')],2);
    drift_phase2_3_within = [drift_phase2_3_within_mean,drift_phase2_3'];
    std_drift_phase2_3 = std(drift_phase2_3_within,0,1,'omitnan')/sqrt(8);
    errorbar(mean(drift_phase2_3_within,'omitnan'),std_drift_phase2_3,'LineWidth',1.5,'Color',[0.72941      0.24706      0.24706]);
    legend([{'Baseline'},{'Novel env'}])
    legend boxoff
    if c==2
        ylim([0.45 1])
    else
        ylim([0.35 0.85])
    end
    ylabel(corr_var_label{c})
    xlim([0 7])
    xticks(1:6)
    xticklabels(0:2:10)
    xlabel('Elapsed time (days)')
    box off
    axis square

   % check sig bw manipulation:
   bs_m = drift_phase1_within';
   enr_m = drift_phase2_3_within';
   [p,tbl,~] = friedman([bs_m(:) enr_m(:)],6);
   p_val_bw_manipulation = round(p,3);
   chi_bw_manipulation = round(tbl{2,5},3);

   disp(['Comparing drift in' corr_var_label{c} 'of baseline vs exposure to novel environments when averaged number of laps is consistent across days'])
   disp(['Using the nonparametric Friedman test: p-value = ' num2str(p_val_bw_manipulation) ' Chi = ' num2str(chi_bw_manipulation)])
end

%% Figure S3B-C, Left - Comparing spatial tuning and activity rates across all pairs of sessions regardless of elapsed time in subsampled data (same number of laps across days):

% Upper panel - Tuning curve correlations
% Lower panel - Ensemble rate correlations 

mean_across = mean(mean_all_sessions_lap_control,2,'omitnan');
clc
for c = 1:2
    mean_phase1 = squeeze(mean_across(1,:,c,:));
    mean_phase2_3 = mean([squeeze(mean_across(2,:,c,:)) ,squeeze(mean_across(3,:,c,:))],2,'omitnan');
    bs_phases = [mean_phase1,mean_phase2_3];
    std_across_mice = std(bs_phases,'omitnan')./sqrt(8);
    figure
    hold on
    plot(bs_phases','color',[0.6 0.6 0.6])
    errorbar(mean(bs_phases,'omitnan'),std_across_mice,'color','k','LineWidth',1.5)
    if c==2
        ylim([0.45 1])
    else
        ylim([0.35 0.85])
    end
    xlim([0.5 2.5])
    ylabel(corr_var_label{c})
    set(gca,'xtick',1:2,'xticklabel',[{'Baseline'},{'Novel env'}])
    box off
    axis square

    % check sig wilcoxson:
    [p_val,~,~] = ranksum(bs_phases(:,1),bs_phases(:,2));
    p_val_wilcxon = round(p_val,3);

    disp(['Comparing ' corr_var_label{c} ' across all pairs of sessions regardless of elapsed time of baseline vs exposure to novel environments when averaged number of laps is consistent across days'])
    disp(['Using the nonparametric Wilcoxon signed-rank test: p-value = ' num2str(p_val_wilcxon)])

end

