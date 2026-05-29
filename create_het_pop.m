% Andy Edwards (Nov 20 2025)

clear; clc; close all

restoredefaultpath;
home = pwd;
addpath(genpath('Baseline models'));
addpath(genpath('Create Population'));
addpath(genpath('Optimization'));
addpath(genpath('Shared'));

% User selections
nTrials = 3000;  % Number of cells in final population
base_paramFile = @Riz2014_init_parameters_INa_low; % Baseline parameter file
Experiments = {'Combined'}; % Accepted as a vector. Possible components: 'Combined','Na','Tail', indicating different optimizations. 'Combined' is published version.
frac_high = [0.370]; % Also accepted as vector, 0.370 is the published optimal population 

for Experiment = 1:length(Experiments)
    %% Specify the 5 best results from the optimization for generating populations
    optimization = Experiments{Experiment}; % 'Combined' or 'Na'
    load([optimization,'_refined.mat']);
    if strcmp(optimization, 'Combined')
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
    
    for i = 1:length(frac_high)
        frac_low = round(1-frac_high(i),3);   

        trials_low = round(nTrials*frac_low);
        trials_high = nTrials - trials_low;

        %Stdev values optimized by DE + experimentally known stdev for other parameters
        if strcmp(optimization,'Combined')
            modParam_stdev_noNa = [0.71 0.90 0.37 0.80 lowestParams(i).g_CaL lowestParams(i).g_CaPQ lowestParams(i).g_CaT 0.94]; %[V_GK_max K_GK g_Kv g_BK g_CaL g_CaPQ g_CaT g_KATP_hat]
        elseif strcmp(optimization,'Na')
            modParam_stdev_noNa = [0.71 0.90 0.37 0.80 0.22 0.59 0.40 0.94]; %[V_GK_max K_GK g_Kv g_BK g_CaL g_CaPQ g_CaT g_KATP_hat]
        end

        % Optimized INa parameters for INa_high and INa_low
        INa_high_modParam_stdev = [lowestParams(i).V_hNa 0.15 lowestParams(i).g_Na lowestParams(i).V_mNa]; %[V_hNa n_hNa g_Na V_mNa]
        INa_low_modParam_stdev = [lowestParams(i).V_hNa_low 0.15 lowestParams(i).g_Na_low lowestParams(i).V_mNa_low]; %[V_hNa_low n_hNa_low g_Na_low V_mNa_low]

        disp(['Rank ',num2str(i),' optimal population construction initiated...'])

        %% Define Models to Analyze
        model = cell(0);
        model = [model, {{base_paramFile, @Riz2014_init_states_INa_low, @Riz2014_rhs_INa_low, 'Beta cell'}}];

        %% Define non INa parameter variation
        modParam_names_noNa = {'V_GK_max', 'K_GK', 'g_Kv', 'g_BK', 'g_CaL', 'g_CaPQ', 'g_CaT', 'g_KATP_hat'};
        PQ_ind = find(strcmp(modParam_names_noNa,'g_CaPQ'));
        nModParams = numel(modParam_names_noNa);
        modParam_scaling_noNa = getScalingFactors(modParam_stdev_noNa, nModParams, nTrials, PQ_ind);

        %% Define INa_high and INa_low parameter variation
        INa_high_modParam_names = {'V_hNa', 'n_hNa', 'g_Na', 'V_mNa'};
        INa_high_nModParams = numel(INa_high_modParam_names);
        INa_high_modParam_scaling = getScalingFactors_INa_inact(INa_high_modParam_stdev, INa_high_modParam_names, trials_high);
        INa_high_modParam_scaling = [INa_high_modParam_scaling; zeros(trials_low,INa_high_nModParams)];

        INa_low_modParam_names = {'V_hNa_low','n_hNa_low', 'g_Na_low', 'V_mNa_low'};
        INa_low_nModParams = numel(INa_low_modParam_names);
        INa_low_modParam_scaling = getScalingFactors_INa_inact(INa_low_modParam_stdev, INa_low_modParam_names, trials_low);
        INa_low_modParam_scaling = [zeros(trials_high,INa_low_nModParams);INa_low_modParam_scaling];

        modParam_scaling = [modParam_scaling_noNa,INa_high_modParam_scaling,INa_low_modParam_scaling];
        modParam_names = [modParam_names_noNa,INa_high_modParam_names, INa_low_modParam_names];
        modParam_stdev = [modParam_stdev_noNa,INa_high_modParam_stdev,INa_low_modParam_stdev];
        nModParams = nModParams+INa_high_nModParams+INa_low_nModParams;

        model_paramFun = model{1}{1};
        model_stateFun = model{1}{2};
        model_rhsFun   = model{1}{3};
        model_name     = model{1}{4};

        %% Create the modified parameter matrix from the scaling factor matrix
        [param_vals_baseline, param_names] = model_paramFun();
        [param_vals_scaled, modParam_baseline, modParam_vals, modParam_inds] = modifyParams(...
            param_vals_baseline, param_names, modParam_scaling,...
            modParam_names, nTrials);

        %% Retrieve baseline state values
        [state_vals_baseline, state_names] = model_stateFun(); 
        %state_vals_baseline has all the inital conditions for 15 isletparams

        %% Initialize output variables
        EE = zeros(nTrials, 1); TE = zeros(nTrials, 1); 
        peak_INa = zeros(nTrials, 1); v_half_act = zeros(nTrials, 1); vha_r2 = zeros(nTrials, 1); peak_ICa = zeros(nTrials, 1); late_ICa = zeros(nTrials, 1);
        v_half = zeros(nTrials, 1); vh_r2 = zeros(nTrials, 1);

        %% Compute output metrics of the population
        parfor iTrial = 1:nTrials 
            state_vals_output = state_vals_baseline';
            [EE(iTrial, 1), ~, TE(iTrial,1)] = Exocytosis(param_vals_scaled(iTrial,:), param_names, state_vals_output, state_names, model_rhsFun);
            [~, v_half_act(iTrial, 1),vha_f(iTrial).mdl, vha_r2(iTrial,1), ~,~] = PeakCurrents(param_vals_scaled(iTrial,:), param_names, state_vals_output, state_names, model_rhsFun);
            [v_half(iTrial, 1), vh_f(iTrial).mdl, vh_r2(iTrial, 1)] = INa_inact(param_vals_scaled(iTrial,:), param_names, state_vals_output, state_names, model_rhsFun);
        end

        %% Save the parameter distributions and output metrics
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

        % Define and create the output folder structure
        outfolder = [home,filesep,'Create Population',filesep,optimization,filesep,sprintf('%.0f',nTrials),filesep,sprintf('%.3f%s',frac_high(i))];
        if ~isfolder(outfolder)
            mkdir(outfolder)
        end
        if ~isempty(sub)
            outfolder = [outfolder,filesep,sub];
            mkdir(outfolder);
        end

        % Save the output metrics
        scaling_name = sprintf('modParam_scaling_%.0f_%.3f.mat',nTrials,frac_high(i)); save([outfolder,filesep,scaling_name], 'modParam_scaling');
        stdev_name = sprintf('modParam_stdev_%.0f_%.3f.mat', nTrials, frac_high(i)); save([outfolder,filesep,stdev_name], 'modParam_stdev');
        param_name = sprintf('modParam_names_%.0f_%.3f.mat', nTrials, frac_high(i)); save([outfolder,filesep,param_name], 'modParam_names');
        peak_INa_name = sprintf('peak_INa_%.0f_%.3f.mat', nTrials, frac_high(i)); save([outfolder,filesep,peak_INa_name], 'peak_INa');
        v_half_act_name = sprintf('v_half_act_%.0f_%.3f.mat', nTrials, frac_high(i)); save([outfolder,filesep,v_half_act_name], 'v_half_act');
        peak_ICa_name = sprintf('peak_ICa_%.0f_%.3f.mat',nTrials, frac_high(i)); save([outfolder,filesep,peak_ICa_name], 'peak_ICa');
        late_ICa_name = sprintf('late_ICa_%.0f_%.3f.mat', nTrials, frac_high(i)); save([outfolder,filesep,late_ICa_name], 'late_ICa');
        EE_name = sprintf('EE_%.0f_%.3f.mat', nTrials, frac_high(i)); save([outfolder,filesep,EE_name], 'EE');
        TE_name = sprintf('TE_%.0f_%.3f.mat', nTrials, frac_high(i)); save([outfolder,filesep,TE_name], 'TE');
        v_half_name = sprintf('v_half_%.0f_%.3f.mat', nTrials, frac_high(i)); save([outfolder,filesep,v_half_name], 'v_half');
        frac_high_name = sprintf('frac_high_%.0f_%.3f.mat', nTrials, frac_high(i)); save([outfolder,filesep,frac_high_name], 'frac_high');

        disp(['Rank ',num2str(i),sub,' population generated successfully.'])

        % Plotting parameter distributions
        % These are the parameters that will be plotted
        param_selection = ["V_GK_max" "K_GK" "g_Kv" "g_BK" "g_CaL" "g_CaPQ" "g_CaT" "g_KATP_hat" "V_hNa" "n_hNa" "g_Na" "V_mNa" "V_hNa_low" "n_hNa_low" "g_Na_low" "V_mNa_low"]; 

        % Reducing the modParam_names, modParam_inds, plot_names, and plot_units variables to those that have been selected.
        [modParam_names,modParam_inds,plot_names,plot_units] = select_plot_params_verbose(param_selection, modParam_names, modParam_inds);
        plotScalingFactors(modParam_scaling, plot_names, nModParams, outfolder);
    end
    clearvars -except home nTrials base_paramFile Experiments frac_high
end
