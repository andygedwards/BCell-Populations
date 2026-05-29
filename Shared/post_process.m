function [] = post_process(loaddir, resultsdir, plotting)

        if ~isfolder(resultsdir)
            mkdir(resultsdir);
            mkdir(strcat(resultsdir,'Figures'));
        end
        %% Initialization of analysis variables
        class_low = cell(nTrials,1);
        class_high = cell(nTrials,1); 
        outs_low = cell(nTrials,1);
        outs_high = cell(nTrials,1);
        param_low = cell(nTrials,1);
        param_high = cell(nTrials,1);
        low_count = 0;
        high_count = 0;
        
        for i = 1:nTrials 
            disp(['Cell number: ',num2str(i)])
            low_file = strcat(celldir,num2str(i),'_low.mat');
            low = load(low_file);
            high_file = strcat(celldir,num2str(i),'_high.mat');
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
       
            if ~isempty(plotting.number)
                if plotting.number == i
                    Title = '';
                    plotTS(low.T,high.T,V_low,V_high,pks_low,pks_high,trs_low,trs_high,i,Title,plot_peaks_troughs)
                end
            elseif ~isempty(plotting.pheno) 
                for j = 1: length(plotting.pheno)
                    if strcmp(plotting.pheno{j},class_low{low_count}) && strcmp(plotting.pheno{j},class_high{high_count})
                        Title = [plotting.pheno{j}, ' in low and high glucose'];
                    elseif strcmp(plotting.pheno{j},class_low{low_count}) 
                        Title = [plotting.pheno{j}, ' in low glucose'];
                    elseif strcmp(plotting.pheno{j},class_high{high_count})
                        Title = [plotting.pheno{j}, ' in high glucose'];
                    else
                        continue
                    end
                    plotTS(low.T,high.T,V_low,V_high,pks_low,pks_high,trs_low,trs_high,i,Title,plot_peaks_troughs)
                end
            elseif ~isempty(plotting.transitions)
                for j = 1: length(plotting.transitions)
                    if strcmp(class_low{low_count},'Silent')
                        if strcmp(plotting.transitions{j},'SiSi') && strcmp(class_high{high_count},'Silent')
                            Title = 'Silent to Silent transition';
                            plotTS(low.T,high.T,V_low,V_high,pks_low,pks_high,trs_low,trs_high,i,Title,plot_peaks_troughs)
                        elseif strcmp(plotting.transitions{j},'SiB') && strcmp(class_high{high_count},'Bursting')
                            Title = 'Silent to Bursting transition';
                            plotTS(low.T,high.T,V_low,V_high,pks_low,pks_high,trs_low,trs_high,i,Title,plot_peaks_troughs)
                        elseif strcmp(plotting.transitions{j},'SiSp') && strcmp(class_high{high_count},'Spiking')
                            Title = 'Silent to Spiking transition';
                            plotTS(low.T,high.T,V_low,V_high,pks_low,pks_high,trs_low,trs_high,i,Title,plot_peaks_troughs)
                        elseif strcmp(plotting.transitions{j},'SiD') && strcmp(class_high{high_count},'Depolarized')
                            Title = 'Silent to Depolarized transition';
                            plotTS(low.T,high.T,V_low,V_high,pks_low,pks_high,trs_low,trs_high,i,Title,plot_peaks_troughs)
                        elseif strcmp(plotting.transitions{j},'SiO') && strcmp(class_high{high_count},'Other')
                            Title = 'Silent to Other transition';
                            plotTS(low.T,high.T,V_low,V_high,pks_low,pks_high,trs_low,trs_high,i,Title,plot_peaks_troughs)
                        else 
                            continue
                        end
                    else
                        continue
                    end
                end
            end
        end
       
        save([resultsdir,filesep,'high_glucose.mat'],'class_high','outs_high','param_high');
        save([resultsdir,filesep,'low_glucose.mat'],'class_low','outs_low','param_low');
        
        clearvars -except home frac_high folder resultsdir all_load_paths all_results_paths nTrials plotting plot_peaks_troughs 
        load([folder,filesep,'metadata',filesep,'metadata.mat'])
        save([resultsdir,filesep,'metadata.mat'], 'modParam_names','nModParams','nTrials','model','modParam_inds')
end