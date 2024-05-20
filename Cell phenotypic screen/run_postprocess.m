% enter the name of the folder you want to postprocess
%typically of the format outputyyyymmddTxxxxx/
clear all
close all
clc

exp_name = 'Same'; %Enter the experiment you want to analyze.
nTrials = 3000; %Enter the number of runs you want to analyze.
frac_high = 1.00; 
folder = sprintf('%s/%.0f_%.2f',exp_name,nTrials,frac_high);
filename = strcat('data/',folder,'/cells/');

resultsdir = strcat('results/',folder);
if ~isfolder(resultsdir)
    mkdir(resultsdir);
    mkdir(strcat(resultsdir,'/Figures'));
end

class_low = cell(nTrials,1);
class_high = cell(nTrials,1); 
outs_low = cell(nTrials,1);
outs_high = cell(nTrials,1);
param_low = cell(nTrials,1);
param_high = cell(nTrials,1);
low_count = 0;
high_count = 0;

for i = 1:nTrials 
    i
    low_file = strcat(filename,num2str(i),'_low.mat');
    low = load(low_file);
    high_file = strcat(filename,num2str(i),'_high.mat');
    high = load(high_file);
    V_low = real(low.Y(:,15));
    V_high = real(high.Y(:,15));
    
    low_count = low_count + 1;
    th = -40;
    [pks_low,pheight_low] = peakfinder(V_low, 10, th, 1, 1);
    [trs_low,trheight_low] = peakfinder(V_low, 10, th, -1, 1);
    [class_low{low_count},outs_low{low_count},param_low{low_count}] = Andrean(pks_low,trs_low,V_low,low.T,low.Y,low.param);
    
    high_count = high_count + 1;
    th = -40;
    [pks_high,pheight_high] = peakfinder(V_high, 10, th, 1, 1);
    [trs_high,trheight_high] = peakfinder(V_high, 10, th, -1, 1); 
    [class_high{high_count},outs_high{high_count},param_high{high_count}] = Andrean(pks_high,trs_high,V_high,high.T,high.Y,high.param);

%     if strcmp(class_high{high_count},'Depolarized') || strcmp(class_low{low_count},'Depolarized')
%         figure
%         subplot(2,1,1), plot(low.T/1000,V_low), hold on, %plot(low.T(pks_low)/1000, V_low(pks_low),'r*', low.T(trs_low)/1000, V_low(trs_low),'g*')
%         ylim([-100,30])
%         ylabel('V_{m} (mV)')
%         title('2 mM glucose', 'Fontsize',22)
%         subplot(2,1,2), plot(high.T/1000,V_high), hold on, %plot(high.T(pks_high)/1000, V_high(pks_high),'r*', high.T(trs_high)/1000, V_high(trs_high),'g*')
%         ylim([-100,30])
%         ylabel('V_{m} (mV)')
%         xlabel('time (s)')
%         title('20 mM glucose', 'Fontsize',22)
%         pause
%         close all
%     end
end

save(strcat(resultsdir,'/high_glucose.mat'),'class_high','outs_high','param_high');
save(strcat(resultsdir,'/low_glucose.mat'),'class_low','outs_low','param_low');

clearvars -except folder resultsdir
load(strcat('data/',folder,'/metadata/metadata.mat'))
save(strcat(resultsdir,'/metadata.mat'), 'modParam_names','nModParams','nTrials','model','modParam_inds')

function [class, outs, param] = Andrean(pks,troughs,V,T,Y,param)
    outs = cell(6,1);
    if isempty(pks)
        class = 'Silent';
        outs{1} = NaN;
        outs{2} = NaN;
        outs{3} = NaN;
        outs{4} = NaN;
        outs{5} = mean(Y(:,14));
        outs{6} = max(Y(:,14));
    else 
        if length(troughs)>length(pks)
            average = mean(V(pks)-V(troughs(1:length(pks))));
        elseif length(troughs)<length(pks)
            if length(pks)-length(troughs) > 1
                average = mean(V(pks(1:length(troughs)))-V(troughs));
            else
                average = mean(V(pks(1:end-1))-V(troughs));
            end  
        else 
            average = mean(V(pks)-V(troughs));
        end
        if average > 10
            [crosses,dummy] = peakfinder(V, 10, -50, -1, 1);
            if isempty(crosses)
                class = 'Depolarized';
                outs{1} = NaN;
                outs{2} = NaN;
                outs{3} = NaN;
                outs{4} = NaN;
                outs{5} = mean(Y(:,14));
                outs{6} = max(Y(:,14));
            elseif length(pks) < 2
                class = 'Depolarized';
                outs{1} = NaN;
                outs{2} = NaN;
                outs{3} = NaN;
                outs{4} = NaN;
                outs{5} = mean(Y(:,14));
                outs{6} = max(Y(:,14));
            else
                maxint = max(diff(T(pks)));
                if maxint > 2000 && numel(pks) > 25
                    diffs = diff(T(pks));
                    end_pks = pks(find(diffs>2000));
                    start_pks = zeros(length(end_pks)+1,1);
                    end_pks(end+1) = pks(end);
                    for i = 2:length(end_pks)
                        start_pks(i) = pks(find(pks==end_pks(i-1))+1);
                    end
                    start_pks(1) = pks(1);
                    for i = 1:length(start_pks)
                        burst_dur(i) = mean(T(end_pks(i))-T(start_pks(i)));
                        startind = find(pks==start_pks(i));
                        endind = find(pks==end_pks(i));
                        burst_freq(i) = 1000/mean(diff(T(pks(startind:endind))));
                    end
                    num_bursts = length(start_pks);
                    burst_dur = mean(burst_dur);
                    burst_freq = mean(burst_freq);                
                    class = 'Bursting';
                    outs{1} = burst_freq;
                    outs{2} = burst_dur;
                    outs{3} = num_bursts;
                    outs{4} = NaN;
                    outs{5} = mean(Y(:,14));
                    outs{6} = max(Y(:,14));
%                     plot(T,V), hold on, plot(T(start_pks),V(start_pks),'r*',T(end_pks),V(end_pks),'g*')
%                     pause
%                     close all
                else
                    if numel(pks) < 50
                        class = 'Other';
                        outs{1} = NaN;
                        outs{2} = NaN;
                        outs{3} = NaN;
                        outs{4} = NaN;
                        outs{5} = mean(Y(:,14));
                        outs{6} = max(Y(:,14));
                    else
                        spike_freq = 1000/mean(diff(T(pks)));
                        class = 'Spiking';
                        outs{1} = NaN;
                        outs{2} = NaN;
                        outs{3} = NaN;
                        outs{4} = spike_freq;
                        outs{5} = mean(Y(:,14));
                        outs{6} = max(Y(:,14));
                    end
                end
            end
        end
    end
end


