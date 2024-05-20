% Code last modified: Andy (6 March 2024)

clear all; clc; close all

restoredefaultpath;
addpath('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/Run_Genetic_Algorithm/V_clamp');
addpath('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/Run_Genetic_Algorithm/PoM')
addpath('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/SKNM/Riz_2014_INa_low')

frac_high = 1.00;
nTrials = 3000;  %Number of cells in final population
exp_name = "Same";
%%
%Final populations
load(sprintf('%s/%.0f_%.2f/modParam_scaling_%.0f_%.2f%.mat',exp_name,nTrials,frac_high,nTrials,frac_high))
modParam_names = {'V_{GKmax}', 'K_{GK}', 'V_{PFKmax}', 'K_{PFK}', 'h_{PFK}', 'V_{GAPDHmax}', 'g_{Kv}', 'g_{BK}', 'g_{CaL}', 'g_{CaPQ}', 'g_{CaT}', 'g_{KATP}', 'g_{HERG}', 'V_{hNa}', 'n_{hNa}', 'g_{Na}', 'V_{mNa}', 'V_{hNalow}','n_{hNalow}', 'g_{Nalow}', 'V_{mNalow}'};
modParam_names = convertCharsToStrings(modParam_names);
nModParams = 21;

figure(1); set(gcf, 'color', 'w'); hold on;
plotind = 0;
for iModParam = 1:nModParams  
    if iModParam == 13
        continue
    end
    plotind = plotind+1;
    subplot(4, ceil(nModParams/4), plotind);
    iScalingDistribution = modParam_scaling(:,iModParam);
    plotDist = iScalingDistribution>0;
    histogram(iScalingDistribution(plotDist), 'Normalization','Probability', 'FaceColor', 'blue');
    title(sprintf('%s', convertCharsToStrings(modParam_names{1,iModParam})))%, mean(iScalingDistribution(plotDist)), std(iScalingDistribution(plotDist)))
%         '\n \x03BC = %.2f; \x03C3 = %.2f' ...
    if plotind == 1 || plotind == 7 || plotind == 13 || plotind == 19 
        ylabel('Fraction of cells')
    end
    if iModParam > (nModParams-6)
        xlabel('Fraction of mean')
    end
    grid on;
end
% sgtitle('Scaling Factor Distributions', 'Interpreter', 'none','FontSize',20);
sgtitle('Fitted Parameter Variation', 'Interpreter', 'none','FontName','Arial','FontSize',30);
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',14, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
figure(1); set(gcf, 'Units', 'Inches', 'Position', [0, 0, 18, 18], 'PaperUnits', 'Inches', 'PaperSize', [18, 18])
%f = gcf; exportgraphics(f,'distribution.png','Resolution', 300)
