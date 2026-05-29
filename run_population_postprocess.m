
    clear; close all; clc

    restoredefaultpath;
    home = pwd;
    addpath(genpath('Glucose screen'));
    addpath(genpath('Baseline models'));
    addpath(genpath('Create Population'));
    addpath(genpath('Optimization'));
    addpath(genpath('Shared'));

    %% Specify the populations that you have simulated and would like to analyze.
    %% If frac_high is left empty, the processing will be performed on all subpopulations that have been analyzed for nTrials
    nTrials = 3000;
    init_frac_high = [0.37]; % Specify the frac_high folders you'd like to analyze, or leave empty to analyze all folders within input_home (below).
    ka = 0.00004; % Defining the source folder for the electro-metabolic coupling variant you would like to analyze
    Experiment = 'Combined';
    input_home = [home,filesep,'Glucose screen',filesep,'data',filesep,Experiment,filesep,sprintf('%.0f',nTrials)];
    results_home = [home,filesep,'Glucose screen',filesep,'results',filesep,Experiment,filesep,sprintf('%.0f',nTrials)];

    %% PLOTTING representative data.
    %% Choose whether to plot the identfied peaks and troughs for the analysis together with the time-series
    plot_peaks_troughs = 0; % boolean

    %% Choose only one of the following 3 plotting options, or leave all empty to suppress representative plotting.
    %% Note that you must save the plots manually after they are presented.

    %% Specify chosen cell numbers to plot (must be in the range 1:nTrials)
    plotting.number = []; % vector with indices in range(1:nTrials) for specific cells you wold like to plot high and low glucose responses
    %% Specify specific phenotypes to plot
    plotting.pheno = {}; % cell array of strings with possible entries {'Silent', 'Bursting', 'Spiking', 'Depolarized', 'Other'}.
    % Using this will plot both high and low glucose time series for any cell with the specified phenotypes either during high or low glucose. 
    %% Specify specific phenotype transitions to plot
    plotting.transitions = {}; % cell array of strings with possible entries {'SiSi, SiB, SiSp, SiD, SiO'}. 
    % Using this will plot both high and low glucose time series for any cell with the following transition pattern from low to high glucose:
    % SiSi: Silent to Silent
    % SiB: Silent to Bursting
    % SiSp: Silent to Spiking
    % SiD: Silent to Depolarized
    % SiO: Silent to Other
    %% Specify whether to allow the user to correct the classification of the plotted phenotypes
    get_user_input = 1;

    % Define the load and results directory structures and assign the frac_high variable
    % accounting for any duplicate frac_high values among the input folders
    % optimization solutions
    load_dir = dir(input_home);
    if isempty(init_frac_high) % Case: no folders are specified, all will be analyzed
        folder_count = 1;
        for ind = 1:length(load_dir)
            if ~strcmp(load_dir(ind).name(1),'.')
                path = [load_dir(ind).folder,filesep,load_dir(ind).name,filesep,sprintf('%s_%.5f','ka',ka)];
                subdir = dir(path);
                isSubfolder = [subdir.isdir] & ~ismember({subdir.name}, {'.', '..'});
                subfolders = {subdir(isSubfolder).name};
                subfolders = subfolders(~startsWith(subfolders, '.'));
                excludeNames = {'cells', 'metadata'};
                subfolders = subfolders(~ismember(subfolders, excludeNames));
                if ~isempty(subfolders)
                    for subind = 1: length(subfolders)
                        all_load_paths{folder_count} = [path,filesep,subfolders{subind}];
                        all_results_paths{folder_count} = [results_home,filesep,load_dir(ind).name,filesep,sprintf('%s_%.5f','ka',ka),filesep,subfolders{subind}];
                        frac_high(folder_count) = str2num(load_dir(ind).name);   
                        folder_count = folder_count + 1;
                    end
                else
                    all_load_paths(folder_count) = {path};
                    all_results_paths{folder_count} = [results_home,filesep,load_dir(ind).name,filesep,sprintf('%s_%.5f','ka',ka)];
                    frac_high(folder_count) = str2num(load_dir(ind).name); 
                    folder_count = folder_count+1;
                end 
            end
        end
    else % Case: user-specified frac_high determines folders to be analyzed
        folder_count = 1;
        for ind = 1:length(init_frac_high)
            if isfolder([input_home,filesep,sprintf('%.3f',init_frac_high(ind))])
                path = [input_home,filesep,sprintf('%.3f',init_frac_high(ind)),filesep,sprintf('%s_%.5f','ka',ka)];
                subdir = dir(path);
                isSubfolder = [subdir.isdir] & ~ismember({subdir.name}, {'.', '..'});
                subfolders = {subdir(isSubfolder).name};
                subfolders = subfolders(~startsWith(subfolders, '.'));
                excludeNames = {'cells', 'metadata'};
                subfolders = subfolders(~ismember(subfolders, excludeNames));
                if ~isempty(subfolders)
                    for subind = 1: length(subfolders)
                        all_load_paths{folder_count} = [path,filesep,subfolders{subind}];
                        all_results_paths{folder_count} = [results_home,filesep,sprintf('%.3f',init_frac_high(ind)),filesep,sprintf('%s_%.5f','ka',ka),filesep,subfolders{subind}];
                        frac_high(folder_count) = init_frac_high(ind);
                        folder_count = folder_count + 1;
                    end
                else
                    all_load_paths(folder_count) = {path};
                    all_results_paths{folder_count} = [results_home,filesep,sprintf('%.3f',init_frac_high(ind)),filesep,sprintf('%s_%.5f','ka',ka)];
                    frac_high(folder_count) = init_frac_high(ind);
                    folder_count = folder_count+1;
                end 
            end
        end
    end

    for h = 1:length(all_load_paths)
        %% Input directory for the simulated population and output directory for the results summaries and figures
        folder = all_load_paths{h};
        celldir = [folder,filesep,'cells',filesep];
        resultsdir = [all_results_paths{h},filesep];
        if ~isfolder(resultsdir)
            mkdir(resultsdir);
            mkdir(strcat(resultsdir,'Summary'));
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
            low_file = strcat(celldir,num2str(i),'_low.mat');
            low = load(low_file);
            high_file = strcat(celldir,num2str(i),'_high.mat');
            high = load(high_file);
            
            if ~isfield(low,'Y')
                V_low = low.Y_low(:,15);
                Y = low.Y_low;
                T = low.T_low;
                param = low.param;
                clearvars low
                low.T = T; low.Y = Y; low.param = param;
                save(strcat(celldir,num2str(i),'_low.mat'),'Y','T','param')
                disp(['Saved correct fieldnames for cell: ', num2str(i),'_low.mat'])
            else
                V_low = real(low.Y(:,15));
                a_low = real(low.Y(:,12));
            end
            if ~isfield(high,'Y')
                V_high = high.Y_high(:,15);
                Y = high.Y_high;
                T = high.T_high;
                param = high.param;
                clearvars high
                high.T = T; high.Y = Y; high.param = param;
                save(strcat(celldir,num2str(i),'_high.mat'),'Y','T','param')
                disp(['Saved correct fieldnames for cell: ', num2str(i),'_high.mat'])
            else
                V_high = real(high.Y(:,15));
                a_high = real(high.Y(:,12));
            end
            
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
            
            [metab_low{low_count},metab_high{high_count}] = check_metab(V_low,a_low,low.T,V_high,a_high,high.T,i);
       
            if ~isempty(plotting.number)
                if plotting.number == i
                    Title = '';
                    [class_low{low_count},class_high{high_count}] = plotTS(low.T,high.T,V_low,V_high,a_low,a_high,pks_low,pks_high,trs_low,trs_high,i,Title,plot_peaks_troughs,get_user_input);
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
                   [class_low{low_count},class_high{high_count}] = plotTS(low.T,high.T,V_low,V_high,a_low,a_high,pks_low,pks_high,trs_low,trs_high,i,Title,plot_peaks_troughs,get_user_input);
                end
            elseif ~isempty(plotting.transitions)
                for j = 1: length(plotting.transitions)
                    if strcmp(class_low{low_count},'Silent')
                        if strcmp(plotting.transitions{j},'SiSi') && strcmp(class_high{high_count},'Silent')
                            Title = 'Silent to Silent transition';
                            [~,class_high{high_count}] = plotTS(low.T,high.T,V_low,V_high,a_low,a_high,pks_low,pks_high,trs_low,trs_high,i,Title,plot_peaks_troughs,get_user_input);
                        elseif strcmp(plotting.transitions{j},'SiB') && strcmp(class_high{high_count},'Bursting')
                            Title = 'Silent to Bursting transition';
                            [~,class_high{high_count}] = plotTS(low.T,high.T,V_low,V_high,a_low,a_high,pks_low,pks_high,trs_low,trs_high,i,Title,plot_peaks_troughs,get_user_input);
                        elseif strcmp(plotting.transitions{j},'SiSp') && strcmp(class_high{high_count},'Spiking')
                            Title = 'Silent to Spiking transition';
                            [~,class_high{high_count}] = plotTS(low.T,high.T,V_low,V_high,a_low,a_high,pks_low,pks_high,trs_low,trs_high,i,Title,plot_peaks_troughs,get_user_input);
                        elseif strcmp(plotting.transitions{j},'SiD') && strcmp(class_high{high_count},'Depolarized')
                            Title = 'Silent to Depolarized transition';
                            [~,class_high{high_count}] = plotTS(low.T,high.T,V_low,V_high,a_low,a_high,pks_low,pks_high,trs_low,trs_high,i,Title,plot_peaks_troughs,get_user_input);
                        elseif strcmp(plotting.transitions{j},'SiO') && strcmp(class_high{high_count},'Other')
                            Title = 'Silent to Other transition';
                            [~,class_high{high_count}] = plotTS(low.T,high.T,V_low,V_high,a_low,a_high,pks_low,pks_high,trs_low,trs_high,i,Title,plot_peaks_troughs,get_user_input);
                        else 
                            continue
                        end
                    else
                        continue
                    end
                end
            end
        end
       
        save([resultsdir,filesep,'high_glucose.mat'],'class_high','outs_high','param_high','metab_high');
        save([resultsdir,filesep,'low_glucose.mat'],'class_low','outs_low','param_low','metab_low');

        clearvars -except home frac_high folder resultsdir all_load_paths all_results_paths nTrials plotting plot_peaks_troughs get_user_input g_katp_hat
        load([folder,filesep,'metadata',filesep,'metadata.mat'])
        save([resultsdir,filesep,'metadata.mat'], 'modParam_names','nModParams','nTrials','model','modParam_inds')
    end

    function [user_choice_low,user_choice_high] = plotTS(lowT,highT,V_low,V_high,a_low,a_high,pks_low,pks_high,trs_low,trs_high,cell_index,plot_title,plot_peaks_and_troughs,get_user_input)
       try
            F = figure; set(F, 'Position', [100, 100, 1200, 1000]);
            subplot(2,2,1), plot(lowT/1000,V_low), hold on, 
            if plot_peaks_and_troughs == 1
                plot(lowT(pks_low)/1000, V_low(pks_low),'r*', lowT(trs_low)/1000, V_low(trs_low),'g*')
            end
            ylim([-100,30])
            ylabel('V_{m} (mV)')
            title('2 mM glucose', 'Fontsize',22)
            subplot(2,2,3), plot(lowT/1000,a_low), hold on,
            ylim([0,10.0])
            ylabel('\textit{a} (mM)')

            subplot(2,2,2), plot(highT/1000,V_high), hold on, 
            if plot_peaks_and_troughs == 1
                plot(highT(pks_high)/1000, V_high(pks_high),'r*', highT(trs_high)/1000, V_high(trs_high),'g*')
            end
            ylim([-100,30])
            ylabel('V_{m} (mV)')
            xlabel('time (s)')
            title('20 mM glucose', 'Fontsize',22)
            subplot(2,2,4), plot(highT/1000,a_high), hold on,
            ylim([0,10.0])
            ylabel('\textit{a} (mM)')

            sgtitle(F,['Cell number: ', num2str(cell_index),'. ', plot_title], 'Fontsize',24)
            if get_user_input
                options = {'Silent', 'Spiking', 'Bursting', 'Depolarized'};
                [low_indx, tf_low] = listdlg('ListString', options, 'SelectionMode', 'single', ...
                         'PromptString', 'Select the class for low glucose:', 'ListSize', [160, 100]);
                [high_indx, tf_high] = listdlg('ListString', options, 'SelectionMode', 'single', ...
                         'PromptString', 'Select the class for high glucose:', 'ListSize', [160, 100]);
    
                if tf_low  % If the user made a selection
                    user_choice_low = options{low_indx};
                    fprintf('You selected %s', user_choice_low,' for the low glucose class.');
                else
                    user_choice_low = 'Other';
                    disp(['No low glucose selection made. Class will be returned as: ','Other']);
                end
                if tf_high  % If the user made a selection
                    user_choice_high = options{high_indx};
                    fprintf('You selected %s', user_choice_high,' for the high glucose class.');
                else
                    user_choice_high = 'Other';
                    disp(['No high glucose selection made. Class will be returned as: ','Other']);
                end
            else
                pause
                user_choice_low = 'Other';
                user_choice_high = 'Other';
            end
                close all
       catch
            msg = ['Failed while attempting to plot ', 'Cell number ', num2str(cell_index),':', plot_title];
       end
    end
    
    function [metab_low, metab_high] = check_metab(Vl,al,Tl,Vh,ah,Th,cell)
        
        % figure('Units', 'normalized', 'Position', [0.1 0.55 0.4 0.4]), subplot(2,1,1), plot(Tl,Vl); subplot(2,1,2), plot(Tl,al);
        % title('Low')
        % figure('Units', 'normalized', 'Position', [0.55 0.55 0.4 0.4]), subplot(2,1,1), plot(Th,Vh); subplot(2,1,2), plot(Th,ah);
        % title('High')
        
        if max(al(Tl>50000))>2.0
            disp(strcat('cell number', num2str(cell), ' = metabolism active at low glucose'))
            metab_low = 'active';
        else
            metab_low = 'inactive';
        end
        if max(ah(Th>50000))>2.0
            disp(strcat('cell number', num2str(cell), ' = metabolism active at high glucose'))
            metab_high = 'active';
        else
            metab_high = 'inactive';
        end
    end


