
clear all
close all
clc

restoredefaultpath;
addpath('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/Cell phenotypic screen')

home = '/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/Cell phenotypic screen/results/';
%% ENTER RESULTS FOLDER FOR ANALYSIS
exp_name = 'Same';
nTrials = 3000;
frac_high = 0.21;
%% ----------------
dir = sprintf('/%s/%.0f_%.2f',exp_name,nTrials,frac_high);
results_dir = strcat(home,dir);

cd(results_dir)
load('low_glucose.mat')
load('high_glucose.mat')
load('metadata.mat')
[pmat_low,param_names_low,param_inds_low,high_inds,low_inds] = convert_param(param_low,modParam_names,modParam_inds);
[pmat_high,param_names_high,param_inds_high,high_inds,low_inds] = convert_param(param_high,modParam_names,modParam_inds);

plot_names = {'V_{GKmax}', 'K_{GK}', 'V_{PFKmax}', 'K_{PFK}', 'h_{PFK}', 'V_{GAPDHmax}', 'g_{Kv}', 'g_{BK}', 'g_{CaL}', 'g_{CaPQ}', 'g_{CaT}', 'g_{KATP}', 'g_{hERG}', 'V_{hNa}', 'n_{hNa}', 'g_{Na}', 'V_{mNa}'};
plot_units = {'M/s', 'mM', 'M/s', 'mM', 'dimensionless', 'M/s', 'nS/pF', 'nS/pF','nS/pF','nS/pF','nS/pF','nS/pF','nS/pF','mV','dimensionless','nS/pF','mV' };

low_classes=unique(class_low,'stable')
low_summary=cellfun(@(x) sum(ismember(class_low,x)),low_classes,'un',0)

high_classes=unique(class_low,'stable')
high_summary=cellfun(@(x) sum(ismember(class_high,x)),high_classes,'un',0)

bar(1:2,[cell2mat(low_summary),cell2mat(high_summary)],'stacked');
xticklabels(["2 mM", "20 mM"]);
ylabel("Cells")
legend(low_classes);
% % exportgraphics(gcf,strcat('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/A-BC 2024/Figs/class_summary_',dir(end-3:end),'.pdf'),'ContentType','vector')
% % exportgraphics(gcf,strcat('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/Cell phenotypic screen/results/',dir,'/Figures/class_summary_',dir(end-3:end),'.pdf'),'ContentType','vector')
% 
low_silent = find(strcmp(class_low,'Silent'));
low_bursters = find(strcmp(class_low,'Bursting'));
low_spikers = find(strcmp(class_low,'Spiking'));
low_depolarized = find(strcmp(class_low,'Depolarized'));
low_other = find(strcmp(class_low,'Other'));

high_silent = find(strcmp(class_high,'Silent'));
high_bursters = find(strcmp(class_high,'Bursting'));
high_spikers = find(strcmp(class_high,'Spiking'));
high_depolarized = find(strcmp(class_high,'Depolarized'));
high_other = find(strcmp(class_high,'Other'));

[rsisi,csisi,silent_to_silent] = find(ismember(low_silent,high_silent));
[silent_to_anything] = setdiff(low_silent,high_silent);
[rsib,csib,silent_to_bursting] = find(ismember(low_silent,high_bursters));
[rsis,csis,silent_to_spiking] = find(ismember(low_silent,high_spikers));
[rsid,csid,silent_to_depolarized] = find(ismember(low_silent,high_depolarized));
[rbsp,cbsp,bursting_to_spiking] = find(ismember(low_bursters,high_spikers));
[rbd,cbd,bursting_to_depolarized] = find(ismember(low_bursters,high_depolarized));
[rspb,cspb,spiking_to_bursting] = find(ismember(low_spikers,high_bursters));
[anything_to_silent] = setdiff(high_silent,low_silent);

for i = 1:length(plot_names)
    mean_matrix_low(1,i) = mean(pmat_low(low_silent,param_inds_low(i)));
    mean_matrix_low(2,i) = mean(pmat_low(low_bursters,param_inds_low(i)));
    mean_matrix_low(3,i) = mean(pmat_low(low_spikers,param_inds_low(i)));
    mean_matrix_low(4,i) = mean(pmat_low(low_depolarized,param_inds_low(i)));
    mean_matrix_low(5,i) = mean(pmat_low(low_other,param_inds_low(i)));

    mean_matrix_high(1,i) = mean(pmat_high(high_silent,param_inds_high(i)));
    mean_matrix_high(2,i) = mean(pmat_high(high_bursters,param_inds_high(i)));
    mean_matrix_high(3,i) = mean(pmat_high(high_spikers,param_inds_high(i)));
    mean_matrix_high(4,i) = mean(pmat_high(high_depolarized,param_inds_high(i)));
    mean_matrix_high(5,i) = mean(pmat_high(high_other,param_inds_high(i)));
    
    mean_matrix_trans(1,i) = mean(pmat_high(rsisi,param_inds_high(i)));
    mean_matrix_trans(2,i) = mean(pmat_high(rsib,param_inds_high(i)));
    mean_matrix_trans(3,i) = mean(pmat_high(rsis,param_inds_high(i)));
end

num_sisi = length(rsisi)
num_sa = numel(silent_to_anything)
num_sib = numel(silent_to_bursting)
num_sis = numel(silent_to_spiking)
num_sid = numel(silent_to_depolarized)
num_bsp = numel(bursting_to_spiking)
num_bd = numel(bursting_to_depolarized)
num_spb = numel(spiking_to_bursting)
num_as = numel(anything_to_silent)

bar(1:5,[num_sib,num_sis,num_sid,num_bsp,num_spb]);
legend({"Silent to Bursting", "Silent to Spiking", "Silent to Depolarized", "Bursting to Spiking", "Spiking to Bursting"});
ylabel("Cells")
% exportgraphics(gcf,strcat('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/A-BC 2024/Figs/transitions_',dir(end-3:end),'.pdf'),'ContentType','vector')
% exportgraphics(gcf,strcat('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/Cell phenotypic screen/results/',dir,'/Figures/transitions_',dir(end-3:end),'.pdf'),'ContentType','vector')

high_sisi = find(ismember(rsisi,high_inds));
high_sib = find(ismember(rsib,high_inds));
high_sis = find(ismember(rsis,high_inds));

high_sisi_ratio = length(high_sisi)/length(rsisi);
high_sib_ratio = length(high_sib)/length(rsib);
high_sis_ratio = length(high_sis)/length(rsis);

% % ANOVAs per dependent (low vs. high)
% class_all = [class_low;class_high];
% param_all = [cell2mat(param_low);cell2mat(param_high)];
% [glucose_all{1:3000}] = deal('low');[glucose_all{3001:6000}] = deal('high');glucose_all = glucose_all';

for i = 1:length(param_inds_low)
    aov_low(i).param = param_names_low(i);
    pmat_si_low = pmat_low(low_silent,:);
    pmat_bu_low = pmat_low(low_bursters,:);
    pmat_sp_low = pmat_low(low_spikers,:);
    pmat_low_short = [pmat_si_low;pmat_bu_low;pmat_sp_low];
    y_low = pmat_low_short(:,param_inds_low(i));
    class_low_short = repmat({'Silent'},length(pmat_si_low(:,1)),1);
    class_low_short = [class_low_short;repmat({'Bursting'},length(pmat_bu_low(:,1)),1)];
    class_low_short = [class_low_short;repmat({'Spiking'},length(pmat_sp_low(:,1)),1)];
    tbl_low = table(class_low_short, y_low);
    aov_low(i).anova = anova(tbl_low, 'y_low ~ class_low_short');
    aov_low(i).s = stats(aov_low(i).anova);
    if aov_low(i).s.pValue(1) < 0.05 
        aov_low(i).mult = multcompare(aov_low(i).anova,"class_low_short");
    end

    aov_high(i).param = param_names_high(i);
    pmat_si_high = pmat_high(high_silent,:);
    pmat_bu_high = pmat_high(high_bursters,:);
    pmat_sp_high = pmat_high(high_spikers,:);
    pmat_high_short = [pmat_si_high;pmat_bu_high;pmat_sp_high];
    y_high = pmat_high_short(:,param_inds_high(i));
    class_high_short = repmat({'Silent'},length(pmat_si_high(:,1)),1);
    class_high_short = [class_high_short;repmat({'Bursting'},length(pmat_bu_high(:,1)),1)];
    class_high_short = [class_high_short;repmat({'Spiking'},length(pmat_sp_high(:,1)),1)];
    tbl_high = table(class_high_short, y_high);
    aov_high(i).anova = anova(tbl_high, 'y_high ~ class_high_short');
    aov_high(i).s = stats(aov_high(i).anova);
    if aov_high(i).s.pValue(1) < 0.05 
        aov_high(i).mult = multcompare(aov_high(i).anova,"class_high_short");
    end

    aov_sia(i).param = param_names_high(i);
    pmat_sib = pmat_high(rsib,:);
    pmat_sisp = pmat_high(rsis,:);
    pmat_sisi = pmat_high(rsisi,:);
    pmat_sia = [pmat_sisi;pmat_sib;pmat_sisp];
    param_inds_sia = param_inds_high;
    y_sia = pmat_sia(:,param_inds_sia(i));
    class_sia = repmat({'SiSi'},length(pmat_sisi(:,1)),1);
    class_sia = [class_sia;repmat({'SiB'},length(pmat_sib(:,1)),1)];
    class_sia = [class_sia;repmat({'SiSp'},length(pmat_sisp(:,1)),1)];
    tbl_sia = table(class_sia, y_sia);
    aov_sia(i).anova = anova(tbl_sia, 'y_sia ~ class_sia');
    aov_sia(i).s = stats(aov_sia(i).anova);
    if aov_sia(i).s.pValue(1) < 0.05 
        aov_sia(i).mult = multcompare(aov_sia(i).anova,"class_sia");
    end
end

figure, x0=10; y0=10; width=1000; height=1200;
set(gcf,'position',[x0,y0,width,height])
grouporder = {'Silent','Bursting','Spiking'};
plotviolins(class_low_short,pmat_low_short,param_inds_low,length(param_inds_low),plot_names,plot_units,aov_low,grouporder)
% exportgraphics(gcf,strcat('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/A-BC 2024/Figs/',dir(end-3:end),'_low_violins.pdf'),'ContentType','vector')
% exportgraphics(gcf,strcat('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/Cell phenotypic screen/results/',dir,'/Figures/',dir(end-3:end),'_low_violins.pdf'),'ContentType','vector')
% 
figure, x0=1050; y0=10; width=1000; height=1200;
set(gcf,'position',[x0,y0,width,height])
grouporder = {'Silent','Bursting','Spiking'};
plotviolins(class_high_short,pmat_high_short,param_inds_high,length(param_inds_high),plot_names,plot_units,aov_high,grouporder)
% exportgraphics(gcf,strcat('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/A-BC 2024/Figs/',dir(end-3:end),'_high_violins.pdf'),'ContentType','vector')
% exportgraphics(gcf,strcat('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/Cell phenotypic screen/results/',dir,'/Figures/',dir(end-3:end),'_high_violins.pdf'),'ContentType','vector')
% 
figure, x0=1050; y0=10; width=1000; height=1200;
set(gcf,'position',[x0,y0,width,height])
grouporder = {'SiSi','SiB','SiSp'};
plotviolins(class_sia,pmat_sia,param_inds_sia,length(param_inds_sia),plot_names,plot_units,aov_sia,grouporder)
% exportgraphics(gcf,strcat('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/A-BC 2024/Figs/',dir(end-3:end),'_trans_violins.pdf'),'ContentType','vector')
% exportgraphics(gcf,strcat('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/Cell phenotypic screen/results/',dir,'/Figures/',dir(end-3:end),'_trans_violins.pdf'),'ContentType','vector')

function [p_matrix,p_mat_names,p_mat_inds,high_inds,low_inds] = convert_param(param,param_names,param_inds)
    p_matrix = zeros(length(param),length(param{1}));
    high_count = 0;
    low_count = 0;
    for i = 1:length(param)
%         i
        if param{i}(param_inds(find(strcmp(param_names,'V_hNa'))))==0 && param{i}(param_inds(find(strcmp(param_names,'n_hNa'))))==0 && param{i}(param_inds(find(strcmp(param_names,'g_Na'))))==0 && param{i}(param_inds(find(strcmp(param_names,'V_mNa'))))==0
            param{i}(param_inds(find(strcmp(param_names,'V_hNa'))))=param{i}(param_inds(find(strcmp(param_names,'V_hNa_low'))));
            param{i}(param_inds(find(strcmp(param_names,'V_mNa'))))=param{i}(param_inds(find(strcmp(param_names,'V_mNa_low'))));
            param{i}(param_inds(find(strcmp(param_names,'n_hNa'))))=param{i}(param_inds(find(strcmp(param_names,'n_hNa_low'))));
            param{i}(param_inds(find(strcmp(param_names,'g_Na'))))=param{i}(param_inds(find(strcmp(param_names,'g_Na_low'))));
            param{i}(param_inds(find(strcmp(param_names,'V_hNa_low')))) = 0;
            param{i}(param_inds(find(strcmp(param_names,'V_mNa_low')))) = 0;
            param{i}(param_inds(find(strcmp(param_names,'n_hNa_low')))) = 0;
            param{i}(param_inds(find(strcmp(param_names,'g_Na_low')))) = 0;
            low_count = low_count+1;
            low_inds(low_count) = i;
        else
            high_count = high_count+1;
            high_inds(high_count) = i;
        end
        p_matrix(i,:) = param{i};
    end
    param_inds(find(strcmp(param_names,'V_hNa_low'))) = [];
    param_names(find(strcmp(param_names,'V_hNa_low'))) = [];
    param_inds(find(strcmp(param_names,'V_mNa_low'))) = [];
    param_names(find(strcmp(param_names,'V_mNa_low'))) = [];
    param_inds(find(strcmp(param_names,'n_hNa_low'))) = [];
    param_names(find(strcmp(param_names,'n_hNa_low'))) = [];
    param_inds(find(strcmp(param_names,'g_Na_low'))) = [];
    param_names(find(strcmp(param_names,'g_Na_low'))) = [];
    
    p_mat_inds = param_inds;
    p_mat_names = param_names;
    if low_count == 0
           low_inds = [];
    end
    if high_count == 0
           high_inds = [];
    end
end

function plotviolins(class,param,paraminds,nParams,names,units,aov,grouporder)
    placeind = 1;
    for i = 1:nParams
        if strcmp(names{i},'g_{hERG}')
            continue
        end
        subplot(4,4,placeind), H = violinplot(param(:,paraminds(i)),class,'GroupOrder',grouporder,'ShowMean',true); 
        for h = 1:length(H)
            H(h).MeanPlot.XData = h;
            H(h).MeanPlot.YData = H(h).MeanPlot.YData(1);
            set(H(h).MeanPlot, 'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','k')
        end
        title(names(i));
        ylabel(units{i})
        if isfield(aov(i),'mult')
            if ~isempty(aov(i).mult)
               sigs{i} = find(aov(i).mult.pValue < 0.05);
               if ~isempty(sigs{i})
                   for ii = 1:length(sigs{i})
                       groups{ii} = {string(table2cell(aov(i).mult(sigs{i}(ii),1))),string(table2cell(aov(i).mult(sigs{i}(ii),2)))};
                       p(ii) = aov(i).mult.pValue(sigs{i}(ii));
                   end
                   H = sigstar(groups,p);
               end
            end
        end
        clearvars sigs groups p
        placeind = placeind+1;
    end
end
