% Code last modified: Andy Edwards (November 2025)

clear; close all; clc

restoredefaultpath;
home = pwd;
addpath(genpath('Baseline models'));
addpath(genpath('Create Population'));
addpath(genpath('Optimization'));
addpath(genpath('Shared'));

%% Load the optimization results that were used for generating populations
Optimizations = {'Combined'};   % cell array including 'Combined' and 'Na', or 'None'
nTrials = 3000;
frac_high = [0.37]; % Define the optimized population for simulating (these must first be created by create_het_pop.m) 
ka = 4e-05; % the electro-metabolic coupling parameter for the simulation set (4e-05 is the published value) 
ka_default = 1e-04; % the default electro-metabolic coupling parameter

for Optimization = 1:length(Optimizations)
    %% Specify the 5 best results from the optimization for generating populations
    if ~strcmp(Optimizations{Optimization},'None')
        Experiment = Optimizations{Optimization}; % 'Combined' or 'Na'
        load([Experiment,'_refined.mat']);
        %% The population size you want to load (must be one of the sizes previously specified when running create_het_pop)
        if strcmp(Experiment, 'Combined')
            [sortedRes, sortedIndices] = sort(vertcat(Combined_refined.Total_Error));
            lowestIndices = sortedIndices(1:min(5, length(sortedIndices)));
            lowestParams = horzcat(Combined_refined(lowestIndices).p);
        else
            [sortedRes, sortedIndices] = sort(vertcat(Na_refined.Total_Error));
            lowestIndices = sortedIndices(1:min(5, length(sortedIndices)));
            lowestParams = horzcat(Na_refined(lowestIndices).p);
        end
        if isempty(frac_high)
            frac_high = vertcat(lowestParams.Frachigh); 
        else
            for i = 1:length(frac_high)
                lowestParams(i) = lowestParams(1);
            end
            lowestParams(length(frac_high)+1:end) = [];
        end
    else 
        Experiment = 'Tail';
    end
    
    % Define the parent directory containing the populations
    load_dir = [home,filesep,'Create Population'];
    
    for i = 1:length(frac_high)
        % Load the population and define the directory for saving the simulated
        % glucose responses
        scaling_fname = strcat(sprintf('modParam_scaling_%.0f_%.3f.mat',nTrials,frac_high(i)));
        names_fname = strcat(sprintf('modParam_names_%.0f_%.3f.mat',nTrials,frac_high(i)));
        
        % If there is a duplicate value for frac_high in the ordered list,
        % create a subfolder directory with labels A and B indicating the order
        % based on cost function error. May need to modify this for cases where more than 2
        % optimization results have identical frac_high values to 3 dec places.
        counts = histc(frac_high, unique(frac_high)); 
        duplicateFlags = counts(arrayfun(@(x) find(unique(frac_high) == x, 1), frac_high)) > 1;
        if i == find(duplicateFlags == 1,1,'first')
            sub = 'A';
        elseif duplicateFlags(i) == 1
            sub = 'B';
        else
            sub = '';
        end
        if ~isempty(sub)
            load(strcat(load_dir,filesep,Experiment,filesep,sprintf('%.0f',nTrials),filesep,sprintf('%.3f',frac_high(i)),filesep,sub,filesep,scaling_fname))
            load(strcat(load_dir,filesep,Experiment,filesep,sprintf('%.0f',nTrials),filesep,sprintf('%.3f',frac_high(i)),filesep,sub,filesep,names_fname))
            savedir = [home,filesep,'Glucose screen',filesep,'data',filesep,Experiment,filesep,sprintf('%.0f',nTrials),filesep,sprintf('%.3f',frac_high(i)),filesep,sub,filesep,sprintf('%s_%.5f','ka',ka)];
            disp(['Low and high glucose simulations initiated for frac_high = ', num2str(frac_high(i)),filesep,sub,' population, at nTrials = ',num2str(nTrials)])  
            if strcmp(Experiment,'Tail')
                modParam_scaling = new_modParam_scaling;
            end
        else
            load(strcat(load_dir,filesep,Experiment,filesep,sprintf('%.0f',nTrials),filesep,sprintf('%.3f',frac_high(i)),filesep,scaling_fname))
            load(strcat(load_dir,filesep,Experiment,filesep,sprintf('%.0f',nTrials),filesep,sprintf('%.3f',frac_high(i)),filesep,names_fname))
            savedir = [home,filesep,'Glucose screen',filesep,'data',filesep,Experiment,filesep,sprintf('%.0f',nTrials),filesep,sprintf('%.3f',frac_high(i)),filesep,sprintf('%s_%.5f','ka',ka)];
            disp(['Low and high glucose simulations initiated for frac_high = ', num2str(frac_high(i)),' population, at nTrials = ', num2str(nTrials)])
            if strcmp(Experiment,'Tail')
                modParam_scaling = new_modParam_scaling;
            end
        end
         disp(['Data will be saved to: ', savedir])  
    
        % Create the save directory structure
        if ~ isfolder(savedir)
            mkdir(savedir);
            mkdir(strcat(savedir,filesep,'cells'));
            mkdir(strcat(savedir,filesep,'metadata'));
        elseif ~ isfolder(strcat(savedir,filesep,'cells'))
            mkdir(strcat(savedir,filesep,'cells'));
        elseif ~ isfolder(strcat(savedir,filesep,'metadata'))
            mkdir(strcat(savedir,filesep,'metadata'));
        end
        
        % Define the variable parameters (used for plotting)
        modParam_names = convertCharsToStrings(modParam_names);
        nModParams = length(modParam_names);
        
        %% Define Models to Analyze
        model = cell(0);
        model = [model, {{@Riz2014_init_parameters_INa_low, @Riz2014_init_states_INa_low, @Riz2014_rhs_INa_low, 'Beta cell'}}];
        
        model_paramFun = model{1}{1};
        model_stateFun = model{1}{2};
        model_rhsFun   = model{1}{3};
        model_name     = model{1}{4};
        
        %% Create the modified parameter matrix from the scaling factor matrix
        [param_vals_baseline, param_names] = model_paramFun();
        [param_vals_scaled, modParam_baseline, modParam_vals, modParam_inds] = modifyParams(...
            param_vals_baseline, param_names, modParam_scaling,...
            modParam_names, nTrials);
        
        % Retrieve Baseline State Values
        [state_vals_baseline, state_names] = model_stateFun();
        %state_vals_baseline has all the inital conditions
    
        % Compute glucose responses of the population   
         runTime = tic;
         disp(['Running simulation on Experiment: ',Experiment, ', and rho: ',num2str(frac_high(i))]);
         parfor iTrial = 1:nTrials
            states = state_vals_baseline';
            init_Tstop = 60000;
            param = param_vals_scaled(iTrial,:);
            % disp(['Original ka = ', num2str(param(find(strcmp(param_names, 'k_A'))))])
            param(find(strcmp(param_names, 'G'))) = 2; % 2mM glucose for initialization 
            param(find(strcmp(param_names, 'k_A'))) = param(find(strcmp(param_names, 'k_A')))*(ka/ka_default);
            % disp(['New ka = ', num2str(param(find(strcmp(param_names, 'k_A'))))])
            options = []; % solver options
            
            % initialization run
            [T_init, Y_init] = ode15s(model_rhsFun, [0, init_Tstop], states, options, param);
        
            Tstop = 600000;  % simulation time
        
            % Low glucose run
            [T_low, Y_low] = ode15s(model_rhsFun, [0, Tstop], Y_init(end,:), options, param);
            fname_low = [savedir,filesep,'cells',filesep,num2str(iTrial),'_low.mat'];
            try 
                parsave(fname_low,T_low,Y_low,param)
            catch
                disp(['Error wile attempting to save:', fname_low])
            end
        
            % High glucose run
            param(find(strcmp(param_names, 'G'))) = 20;
            [T_high, Y_high] = ode15s(model_rhsFun, [0, Tstop], Y_init(end,:), options, param);
            fname_high = [savedir,filesep,'cells',filesep,num2str(iTrial),'_high.mat'];
            try
                parsave(fname_high,T_high,Y_high,param)
            catch
                disp(['Error wile attempting to save:', fname_high])
            end
         end
        save(strcat(savedir,filesep,'metadata',filesep,'metadata.mat'),'modParam_names','nModParams','nTrials','model','modParam_inds');
        disp(['Low and high glucose simulations for population: ', num2str(frac_high(i)),' were successful'])
    end
    clearvars -except home nTrials Optimizations frac_high ka ka_default
end
function parsave(fname,T,Y, param)
  save(fname,'T', 'Y','param');
end



