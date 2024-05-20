clear; clc; close all

% Code last modified: Andy Edwards (5th Feb 2024)

M = readtable('beta_glucose180.csv');  %Experimental values from the Camunas-Soler (2019) paper https://github.com/jcamunas/patchseq

num = 3000;
frac_high = 0.00;
pop_folder = strcat('Same/',sprintf('%.0f_%.2f', num, frac_high));
cd(strcat('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/Run_Genetic_Algorithm/Create_optimized_population/',pop_folder));

% sprintf('peak_INa_%.2f_%.2f.mat', num, frac_high)


% Experimental Data
Late_CaL_current = M.LateCa2_Current;
early_ca_current = M.EarlyCa2_Current;
peak_INa_current = M.PeakNa_Current;
half_inact_sodium_current = M.HalfInactivationSodiumCurrent_mV;
total_exo = M.TotalExocitosis;
early_exo = M.EarlyExocytosis;

early_ca_current = early_ca_current(~isnan(early_ca_current));
Late_CaL_current = Late_CaL_current(~isnan(Late_CaL_current));Late_CaL_current = Late_CaL_current(Late_CaL_current<=0);
peak_INa_current = peak_INa_current(~isnan(peak_INa_current));
half_inact_sodium_current = half_inact_sodium_current(~isnan(half_inact_sodium_current));
half_inact_sodium_current = half_inact_sodium_current(half_inact_sodium_current<=0);
total_exo = total_exo(~isnan(total_exo)); total_exo = total_exo(total_exo>=0);
early_exo = early_exo(~isnan(early_exo)); early_exo = early_exo(early_exo>=0);

% Simulated population
peak_INa_sim = load(sprintf('peak_INa_%.0f_%.2f.mat', num, frac_high));
early_CaL_sim = load(sprintf('peak_ICa_%.0f_%.2f.mat', num, frac_high));
late_CaL_sim = load(sprintf('late_ICa_%.0f_%.2f.mat', num, frac_high));
v_half_sim = load(sprintf('v_half_%.0f_%.2f.mat', num, frac_high));
ee_sim = load(sprintf('EE_%.0f_%.2f.mat', num, frac_high));
te_sim = load(sprintf('TE_%.0f_%.2f.mat', num, frac_high));
load(sprintf('modParam_stdev_%.0f_%.2f.mat', num, frac_high));
load(sprintf('modParam_names_%.0f_%.2f.mat', num, frac_high));
load(sprintf('modParam_scaling_%.0f_%.2f.mat', num, frac_high));
load(sprintf('frac_high_%.0f_%.2f.mat', num, frac_high))

high_cols = [find(strcmp(modParam_names,'V_mNa')),find(strcmp(modParam_names,'V_hNa')),find(strcmp(modParam_names,'n_hNa')),find(strcmp(modParam_names,'g_Na'))];
high_rows = find(all(~modParam_scaling(:,high_cols) == 0,2));
low_cols = [find(strcmp(modParam_names,'V_mNa_low')),find(strcmp(modParam_names,'V_hNa_low')),find(strcmp(modParam_names,'n_hNa_low')),find(strcmp(modParam_names,'g_Na_low'))];
low_rows = find(all(~modParam_scaling(:,low_cols) == 0,2));

early_CaL_sim = cell2mat(struct2cell(early_CaL_sim));
late_CaL_sim = cell2mat(struct2cell(late_CaL_sim));
peak_INa_sim = cell2mat(struct2cell(peak_INa_sim));
peak_INa_sim_high = peak_INa_sim(high_rows);
peak_INa_sim_low = peak_INa_sim(low_rows);
v_half_sim = cell2mat(struct2cell(v_half_sim));
v_half_sim_high = v_half_sim(high_rows);
v_half_sim_low = v_half_sim(low_rows);
ee_sim = cell2mat(struct2cell(ee_sim));
te_sim = cell2mat(struct2cell(te_sim));

%%
% Normalize data to mean

figure(1); set(gcf, 'color', 'w'); hold on;

subplot(2, 3, 1); hold on
range_earlyCaL = abs(min(early_CaL_sim./mean(early_CaL_sim)) - max(early_CaL_sim./mean(early_CaL_sim)));
bw1 = range_earlyCaL/19;
h1 = histogram(early_ca_current./mean(early_ca_current),'BinWidth',bw1, 'Normalization', 'Probability','facecolor','b', 'facealpha',.4);
hold on; 
h2 = histogram(early_CaL_sim./mean(early_CaL_sim),'BinWidth',bw1,'Normalization', 'Probability', 'facecolor','r','facealpha',.4);
title('Early {\it{I_{Ca}}}'); 
xlabel('Normalized Early {\it{I_{Ca}}}'); ylabel('Frequency')

subplot(2, 3, 2); hold on
range_lateCaL = abs(min(late_CaL_sim./mean(late_CaL_sim)) - max(late_CaL_sim./mean(late_CaL_sim)));
bw2 = range_lateCaL/19;
histogram(Late_CaL_current./mean(Late_CaL_current),'BinWidth',bw2, 'Normalization', 'Probability','facecolor','b', 'facealpha',.4);
hold on; 
histogram(late_CaL_sim./mean(late_CaL_sim),'BinWidth',bw2,'Normalization', 'Probability', 'facecolor','r','facealpha',.4);
title('Late {\it{I_{Ca}}} '); 
xlabel('Normalized Late {\it{I_{Ca}}} '); ylabel('Frequency')


subplot(2, 3, 3); hold on
range_peakINa = abs(min(peak_INa_sim./mean(peak_INa_sim)) - max(peak_INa_sim./mean(peak_INa_sim)));
bw3 = range_peakINa/19;
histogram(peak_INa_current./mean(peak_INa_current),'BinWidth',bw3, 'Normalization', 'Probability','facecolor','b', 'facealpha',.4);
hold on; 
histogram(peak_INa_sim./mean(peak_INa_sim),'BinWidth',bw3,'Normalization', 'Probability', 'facecolor','r','facealpha',.4);
title('{\it{I_{Na}}} peak'); 
xlabel('Normalized {\it{I_{Na}}} peak'); ylabel('Frequency')


subplot(2, 3, 4); hold on
range_vhalf = abs(min(half_inact_sodium_current) - max(half_inact_sodium_current));
bw4 = range_vhalf/19;
histogram(half_inact_sodium_current,'BinWidth',bw4, 'Normalization', 'Probability','facecolor','b', 'facealpha',.4);
hold on; 
histogram(v_half_sim,'BinWidth',bw4,'Normalization', 'Probability', 'facecolor','r','facealpha',.4);
title(' {\it{I_{Na}}} half inactivation'); 
xlabel('{\it{I_{Na}}} half inactivation (mV)'); ylabel('Frequency')

subplot(2, 3, 5); hold on
range_ee = abs(min(early_exo./mean(early_exo)) - max(early_exo./mean(early_exo)));
bw5 = range_ee/19;
histogram(early_exo./mean(early_exo),'BinWidth',bw5, 'Normalization', 'Probability','facecolor','b', 'facealpha',.4);
hold on; 
histogram(ee_sim./mean(ee_sim),'BinWidth',bw5,'Normalization', 'Probability', 'facecolor','r','facealpha',.4);
title('Early Exocytosis'); 
xlabel('Normalized Early Exocytosis'); ylabel('Frequency')

subplot(2, 3, 6); hold on
range_te = abs(min(total_exo./mean(total_exo)) - max(total_exo./mean(total_exo)));
bw6 = range_te/19;
histogram(total_exo./mean(total_exo),'BinWidth',bw6, 'Normalization', 'Probability','facecolor','b', 'facealpha',.4);
hold on; 
histogram(te_sim./mean(te_sim),'BinWidth',bw6,'Normalization', 'Probability', 'facecolor','r','facealpha',.4);
title('Total Exocytosis'); 
xlabel('Normalized Total Exocytosis'); ylabel('Frequency')
legend('Experimental Distribution', 'Simulated Distribution','NumColumns',2)

%%
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',13.5, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
figure(1); set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 8], 'PaperUnits', 'Inches', 'PaperSize', [10, 8])
%f = gcf; exportgraphics(f,'report_fig_Roshni.png','Resolution', 300)

figure(2); set(gcf, 'color', 'w'); hold on;

subplot(1, 2, 1); hold on
range_peakINa = abs(min(peak_INa_sim./mean(peak_INa_sim)) - max(peak_INa_sim./mean(peak_INa_sim)));
bw3 = range_peakINa/19;
histogram(peak_INa_current./mean(peak_INa_current),'BinWidth',bw3, 'Normalization', 'Probability','facecolor','b', 'facealpha',.4);
hold on; 
% histogram(peak_INa_sim./mean(peak_INa_sim),'BinWidth',bw3,'Normalization', 'Probability', 'facecolor','r','facealpha',.4);
histogram(peak_INa_sim_low./mean(peak_INa_sim),'BinWidth',bw3,'Normalization', 'Probability', 'facecolor','g','facealpha',.4);
histogram(peak_INa_sim_high./mean(peak_INa_sim),'BinWidth',bw3,'Normalization', 'Probability', 'facecolor','m','facealpha',.4);
title('{\it{I_{Na}}} peak'); 
xlabel('Normalized {\it{I_{Na}}} peak'); ylabel('Frequency')


subplot(1, 2, 2); hold on
range_vhalf = abs(min(half_inact_sodium_current) - max(half_inact_sodium_current));
bw4 = range_vhalf/19;
histogram(half_inact_sodium_current,'BinWidth',bw4, 'Normalization', 'pdf','facecolor','b', 'facealpha',.4);
hold on; 
% histogram(v_half_sim,'BinWidth',bw4,'Normalization', 'pdf', 'facecolor','r','facealpha',.4);
histogram(v_half_sim_low,'BinWidth',bw4,'Normalization', 'pdf', 'facecolor','g','facealpha',.4);
histogram(v_half_sim_high,'BinWidth',bw4,'Normalization', 'pdf', 'facecolor','m','facealpha',.4);
title(' {\it{I_{Na}}} half inactivation'); 
xlabel('{\it{I_{Na}}} half inactivation (mV)'); ylabel('Frequency')

set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',13.5, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
figure(2); set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 8], 'PaperUnits', 'Inches', 'PaperSize', [10, 8])

cd('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/Run_Genetic_Algorithm/Create_optimized_population/');

