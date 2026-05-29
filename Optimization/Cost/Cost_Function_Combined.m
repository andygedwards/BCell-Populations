function Total_Error = Cost_Function_Combined(nTrials, P)

    % Code last modified: Andy (October 2024)

    home = pwd;
    addpath(genpath(home))
    exp = readtable(fullfile(home,'Input Data','Cell','beta_glucose180.csv'));  % Load experimental values
    early_ICa_exp = exp.EarlyCa2_Current(~isnan(exp.EarlyCa2_Current));
    late_ICa_exp = exp.LateCa2_Current(~isnan(exp.LateCa2_Current));
    EE_exp = exp.EarlyExocytosis(~isnan(exp.EarlyExocytosis));
    TE_exp = exp.TotalExocitosis(~isnan(exp.TotalExocitosis));
    peak_INa_current = exp.PeakNa_Current(~isnan(exp.PeakNa_Current));
    half_inact_sodium_current = exp.HalfInactivationSodiumCurrent_mV(~isnan(exp.HalfInactivationSodiumCurrent_mV));
    half_inact_sodium_current = half_inact_sodium_current(half_inact_sodium_current<=0);

    % Fraction of cells expressing high voltage INa
    frac_high = round(P.Frachigh,2);  
    frac_low = round(1-frac_high,2);  

    trials_low = round(nTrials*frac_low);
    trials_high = round(nTrials-trials_low);

    %% Non-optimized parameter variation
    st_dev_vector(1) = 0.71; % V_GK_max stdev (direct from experiments)
    st_dev_vector(2) = 0.90; % K_GK stdev (direct from experiments)
    st_dev_vector(3) = 0.37; % g_Kv stdev (direct from experiments)
    st_dev_vector(4) = 0.80; % g_BK stdev (direct from experiments)
    st_dev_vector(8) = 0.94; % g_KATP_hat stdev (direct from experiments)
    st_dev_vector(10) = 0.15; % n_hNa stdev (direct from experiments)
    st_dev_vector(14) = 0.15; % n_hNa_low stdev (direct from experiments)

    %% Parameter variations passed by optimizer
    st_dev_vector(5) = P.g_CaL;
    st_dev_vector(6) = P.g_CaPQ;
    st_dev_vector(7) = P.g_CaT;
    st_dev_vector(9) = P.V_hNa;
    st_dev_vector(11) = P.g_Na;
    st_dev_vector(12) = P.V_mNa;
    st_dev_vector(13) = P.V_hNa_low;
    st_dev_vector(15) = P.g_Na_low;
    st_dev_vector(16) = P.V_mNa_low;

    modParam_stdev_noNa = st_dev_vector(1:8); % standard deviations (normalized to the mean) of the non INa parameters
    INa_high_modParam_stdev = st_dev_vector(9:12); % standard deviations (normalized to the mean) of the INa_high parameters
    INa_low_modParam_stdev = st_dev_vector(13:16); % standard deviations (normalized to the mean) of the INa_low parameters
    
    %% Define Models to Analyze
    model = cell(0);
    model = [model, {{@Riz2014_init_parameters_INa_low, @Riz2014_init_states_INa_low, @Riz2014_rhs_INa_low, 'Beta cell'}}];

    %% Define non INa parameter variation
    modParam_names_noNa = {'V_GK_max', 'K_GK', 'g_Kv', 'g_BK', 'g_CaL', 'g_CaPQ', 'g_CaT', 'g_KATP_hat'};
    PQ_ind = find(strcmp(modParam_names_noNa,'g_CaPQ'));
    nModParams = numel(modParam_names_noNa);
    modParam_scaling_noNa = getScalingFactors(modParam_stdev_noNa, nModParams, nTrials, PQ_ind);
    
    %% Define INa_high and INa_low parameter variation
    INa_high_modParam_names = {'V_hNa', 'n_hNa', 'g_Na', 'V_mNa'};
    INa_high_nModParams = numel(INa_high_modParam_names);
    INa_high_modParam_scaling = getScalingFactors_INa_inact(INa_high_modParam_stdev, INa_high_modParam_names, trials_high);
    INa_high_modParam_scaling = [INa_high_modParam_scaling;zeros(trials_low,INa_high_nModParams)];
    
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
    
    %%
    [param_vals_baseline, param_names] = model_paramFun();
    [param_vals_scaled, modParam_baseline, modParam_vals] = modifyParams(...
        param_vals_baseline, param_names, modParam_scaling,...
        modParam_names, nTrials);
    
    % Retrieve Baseline State Values
    [state_vals_baseline, state_names] = model_stateFun();
    %state_vals_baseline has al l the inital conditions for 15 isletparams
    nStates = numel(state_vals_baseline);
    EE = zeros(nTrials, 1); TE = zeros(nTrials, 1); 
    peak_INa = zeros(nTrials, 1); v_half_act = zeros(nTrials, 1); vha_r2 = zeros(nTrials, 1); peak_ICa = zeros(nTrials, 1); late_ICa = zeros(nTrials, 1);
    v_half = zeros(nTrials, 1); vh_r2 = zeros(nTrials, 1);
    
    % Compute the output metrics for the population
    try
        parfor iTrial = 1:nTrials 
            state_vals_output = state_vals_baseline';
            
            %[EE, LE, TE, v_half, vh_r2, peak_INa, v_half_act, vha_r2, peak_ICa, late_ICa, ramp_IK, Peak_IK_minus20, Peak_IK_max, v_half_act_IK, n_act_IK] = run_v_clamp(param_vals_scaled(iTrial,:), param_names,...
                 %state_vals_output, state_names, model_rhsFun);
            [EE(iTrial, 1), ~, TE(iTrial,1)] = Exocytosis(param_vals_scaled(iTrial,:), param_names, state_vals_output, state_names, model_rhsFun);
            [peak_INa(iTrial, 1), v_half_act(iTrial, 1), vha_r2(iTrial, 1), early_ICa(iTrial, 1), late_ICa(iTrial, 1)] = PeakCurrents(param_vals_scaled(iTrial,:), param_names, state_vals_output, state_names, model_rhsFun);
            [v_half(iTrial, 1), vh_r2(iTrial, 1)] = INa_inact(param_vals_scaled(iTrial,:), param_names, state_vals_output, state_names, model_rhsFun);
        end
        close all
    
        %%  Early exocytosis (EE) error 
        EE_scaling = 0.5;
        EE_raw_error = pearson_distribution(EE_exp, EE, 'Early exocytosis');
        EE_error = EE_scaling.*EE_raw_error;
        
        %% Total Exocytosis (TE) error
        TE_scaling = 0.5;
        TE_raw_error = pearson_distribution(TE_exp, TE, 'Total exocytosis');
        TE_error = TE_scaling.*TE_raw_error;
        
        %% Early ICa error
        early_ICa_scaling = 1;
        early_ICa_raw_error = pearson_distribution(early_ICa_exp, early_ICa, 'Early I_{Ca}');
        early_ICa_error = early_ICa_scaling.*early_ICa_raw_error;
        
        %% Late ICa error
        late_ICa_scaling = 1;
        late_ICa_raw_error = pearson_distribution(late_ICa_exp, late_ICa, 'Late I_{Ca}');
        late_ICa_error = late_ICa_scaling.*late_ICa_raw_error;
        
        %% peak INa error
        peak_INa_scaling = 1;
        peak_INa_raw_error = pearson_distribution(peak_INa_current, peak_INa, 'Peak $I_{Na}$');
        peak_INa_error = peak_INa_scaling.*peak_INa_raw_error;
        
        %% v_half INa error
        vhalf_scaling = 3;
        vhalf_raw_error = double_pearson_distribution(half_inact_sodium_current, v_half, '$I_{Na}$ half inactivation');
        vhalf_error = vhalf_scaling.*vhalf_raw_error;
    
        %% Total error
        Total_Error = EE_error + TE_error + early_ICa_error + late_ICa_error + peak_INa_error + vhalf_error;
    
    catch
        disp('Error encountered');
        Total_Error = 5;
    end

function error = pearson_distribution(exp_data, sim_data, data_name)
        % Normalize the experimental and simulated data by their means
        normalized_exp_data = exp_data ./ mean(exp_data);
        normalized_sim_data = sim_data ./ mean(sim_data);
    
        % Calculate moments for experimental data
        exp_mu = mean(normalized_exp_data);
        exp_sigma = std(normalized_exp_data);
        exp_skewness = skewness(normalized_exp_data);
        exp_kurtosis = kurtosis(normalized_exp_data);
    
        % Calculate moments for simulated data
        sim_mu = mean(normalized_sim_data);
        sim_sigma = std(normalized_sim_data);
        sim_skewness = skewness(normalized_sim_data);
        sim_kurtosis = kurtosis(normalized_sim_data);
    
        % Define common x-values to evaluate the PDFs
        x_values = linspace(min(cellfun(@min,{normalized_exp_data, normalized_sim_data})), max(cellfun(@max,{normalized_exp_data, normalized_sim_data})), 1000);
    
        % PDF normalization
        bin_width = (max(normalized_exp_data) - min(normalized_exp_data)) / 15;  % Define a reasonable bin width

        % Generate Pearson distribution for experimental data
        exp_pdf_values = pearson_pdf(x_values, [exp_mu, exp_sigma, exp_skewness, exp_kurtosis]);
        exp_pdf_values = exp_pdf_values * bin_width;
    
        % Generate Pearson distribution for simulated data
        sim_pdf_values = pearson_pdf(x_values, [sim_mu, sim_sigma, sim_skewness, sim_kurtosis]);
        sim_pdf_values = sim_pdf_values * bin_width;

        % Plot normalized histogram of the experimental data
        bin_width = (max(normalized_exp_data) - min(normalized_exp_data)) / 15;  % Define a reasonable bin width
        figure;
        histogram(normalized_exp_data, 'BinWidth', bin_width, 'Normalization', 'Probability', 'FaceColor', 'b', 'FaceAlpha', 0.2);
        hold on;

        % Plot Pearson PDF for the experimental data
        plot(x_values, exp_pdf_values, 'b-', 'LineWidth', 2); % Experimental data Pearson PDF

        % Plot Pearson PDF for the simulated data
        plot(x_values, sim_pdf_values, 'r-', 'LineWidth', 2); % Simulated data Pearson PDF

        % Add labels and title
        title(['Normalized ', data_name, ' - Histogram and Pearson Distribution']);
        xlabel(['Normalized ', data_name]);
        ylabel('Probability');

        % Display legend
        legend('Experimental Data (Histogram)', 'Experimental Data (Pearson)', 'Simulated Data (Pearson)', 'Location', 'Best');
        
        % Compute the RMS error between the experimental and simulated Pearson PDFs
        error = sqrt(mean((sim_pdf_values - exp_pdf_values).^2));
        
        % Display the RMS error
        disp([data_name, 'RMS Error: ', num2str(error)]);
  
    end

    function error = double_pearson_distribution(exp_data, sim_data, data_name)
        % Normalize the experimental and simulated data by their means
        normalized_exp_data = exp_data ./ mean(exp_data);
        normalized_sim_data = sim_data ./ mean(sim_data);

        bin_width = (max(normalized_exp_data) - min(normalized_exp_data)) / 15;  % Define a reasonable bin width
    
        % Fit two Pearson distributions to experimental data
        [exp_weights, exp_params1, exp_params2] = fit_double_pearson(normalized_exp_data);
    
        % Fit two Pearson distributions to simulated data
        [sim_weights, sim_params1, sim_params2] = fit_double_pearson(normalized_sim_data);
    
        % Define common x-values to evaluate the PDFs
        x_values = linspace(min(cellfun(@min,{normalized_exp_data, normalized_sim_data})), max(cellfun(@max,{normalized_exp_data, normalized_sim_data})), 1000);
    
        % Evaluate Pearson mixture PDF for experimental data
        exp_pdf1_values = exp_weights(1) * pearson_pdf(x_values, exp_params1);
        exp_pdf1_values = exp_pdf1_values * bin_width;
        exp_pdf2_values = exp_weights(2) * pearson_pdf(x_values, exp_params2); 
        exp_pdf2_values = exp_pdf2_values * bin_width;
        exp_pdf_values =  exp_pdf1_values + exp_pdf2_values;           
    
        % Evaluate Pearson mixture PDF for simulated data
        sim_pdf1_values = sim_weights(1) * pearson_pdf(x_values, sim_params1); 
        sim_pdf1_values = sim_pdf1_values * bin_width;
        sim_pdf2_values = sim_weights(2) * pearson_pdf(x_values, sim_params2); 
        sim_pdf2_values = sim_pdf2_values * bin_width;
        sim_pdf_values =  sim_pdf1_values + sim_pdf2_values; 
    
        % Plot the normalized histogram of the experimental data
        figure;
        histogram(normalized_exp_data, 'BinWidth', bin_width, 'Normalization', 'Probability', 'FaceColor', 'b', 'FaceAlpha', 0.2);
        hold on;

        % Plot the Pearson mixture PDF for the experimental data
        plot(x_values, exp_pdf1_values, 'b-', 'LineWidth', 2); % 1st component of experimetnal PDF
        plot(x_values, exp_pdf2_values, 'g-', 'LineWidth', 2); % 2nd component of experimental PDF
        plot(x_values, exp_pdf_values, 'k-', 'LineWidth', 2); % 2nd component of experimental PDF

        % Plot the Pearson mixture PDF for the simulated data
        plot(x_values, sim_pdf1_values, 'y-', 'LineWidth', 2); % 1st component of simulated PDF
        plot(x_values, sim_pdf2_values, 'm-', 'LineWidth', 2); % 2nd component of simulated PDF
        plot(x_values, sim_pdf_values, 'r-', 'LineWidth', 2); % 1st component of simulated PDF

        % Add labels and title
        title(['Normalized ', data_name, ' - Histogram and Pearson Mixture Model']);
        xlabel(['Normalized ', data_name]);
        ylabel('Probability');

        % Display legend
        legend('Experimental Data (Histogram)', 'Experimental component 1 (Pearson Mixture)','Experimental component 2 (Pearson Mixture)', 'Experimental combined (Pearson Mixture)',...
            'Simulated component 1 (Pearson Mixture)', 'Simulated component 2 (Pearson Mixture)', 'Simulated combined (Pearson Mixture)','Location', 'Best');
        
        % Compute the RMS error between the experimental and simulated Pearson mixture PDFs
        error = sqrt(mean((sim_pdf_values - exp_pdf_values).^2));
        
        % Display the RMS error
        disp([data_name, 'RMS Error: ', num2str(error)]);
    
        % Hold off to finish plotting
        hold off;
    end
    
    function [weights, params1, params2] = fit_double_pearson(data)
        % Fit two Pearson distributions to the given data using method of moments
        
        % Initial guesses for [mu, sigma, skewness, kurtosis]
        mean1 = mean(data(data > mean(data)));
        mean2 = mean(data(data < mean(data)));
        std1 = std(data(data > mean(data)));
        std2 = std(data(data < mean(data)));
        skew1 = -abs(skewness(data(data > mean(data))));
        skew2 = -abs(skewness(data(data < mean(data))));
        kurt1 = kurtosis(data(data > mean(data)));
        kurt2 = kurtosis(data(data < mean(data)));
        
        % Ensure initial guesses for kurtosis satisfy the constraint kurt >= skew^2 + 1
        kurt1 = max(kurt1, skew1^2 + 1);
        kurt2 = max(kurt2, skew2^2 + 1);
        
        initial_params1 = [mean1, std1, skew1, kurt1];
        initial_params2 = [mean2, std2, skew2, kurt2];
        
        % Initial guesses for the weights of the two components
        initial_weights = [0.5, 0.5]; % Equal weight for both components
        
        % Combine initial parameters into a single vector for optimization
        initial_guess = [initial_weights, initial_params1, initial_params2];
        
        % Define bounds for the optimization
        lb = [0, 0, -Inf, 0.01, -Inf, 1, -Inf, 0.01, -Inf, 1]; % Lower bounds
        ub = [1, 1, Inf, Inf, 0, Inf, Inf, Inf, 0, Inf]; % Upper bounds
        
        % Objective function to minimize the negative log-likelihood of the data
        objective = @(params) -sum(log(pearson_mixture_pdf(data, params(1:2), params(3:6), params(7:10))));
        
        % Define the nonlinear constraint function
        nonlcon = @(params) nonlinear_constraints(params);
        
        % Optimization options
        options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point');
        
        % Optimize the parameters using fmincon
        optimized_params = fmincon(objective, initial_guess, [], [], [], [], lb, ub, nonlcon, options);
        
        % Extract the optimized weights and parameters
        weights = optimized_params(1:2);
        params1 = optimized_params(3:6); % [mean1, std1, skew1, kurt1]
        params2 = optimized_params(7:10); % [mean2, std2, skew2, kurt2]
    end
    
    function pdf_values = pearson_mixture_pdf(x_values, weights, params1, params2)
        % Calculate the PDF of a mixture of two Pearson distributions
        pdf_values = weights(1) * pearson_pdf(x_values, params1) + weights(2) * pearson_pdf(x_values, params2);
        pdf_values(find(pdf_values<1e-6)) = 1e-6;
    end
    
    function pdf_values = pearson_pdf(x_values, params)
        % Generate the Pearson PDF using the moments: mean, std, skewness, kurtosis
        mu = params(1);
        sigma = params(2);
        skewness_val = params(3);
        kurtosis_val = params(4);
    
        % Generate a Pearson random distribution based on the given parameters
        pearson_data = pearsrnd(mu, sigma, skewness_val, kurtosis_val, 1, 10000);
        
        % Estimate the kernel density of the generated Pearson distribution
        if ~all(isnan(pearson_data)) && ~all(isnan(x_values))
            optimal_bandwidth = std(pearson_data) / 2;
            [pdf_values, ~] = ksdensity(pearson_data, x_values,'Bandwidth', optimal_bandwidth);
        else
            pdf_values = 0;
        end
    end

    % Nonlinear constraint function to ensure kurtosis >= skew^2 + 1
    function [c, ceq] = nonlinear_constraints(params)
        % Extract skewness and kurtosis parameters
        skew1 = params(5);  % Skewness1
        kurt1 = params(6);  % Kurtosis1
        skew2 = params(9);  % Skewness2
        kurt2 = params(10); % Kurtosis2
        
        % Inequality constraints (c <= 0)
        % kurt1 >= skew1^2 + 1 and kurt2 >= skew2^2 + 1
        c = [
            (skew1^2 + 1) - kurt1; % Ensure kurt1 >= skew1^2 + 1
            (skew2^2 + 1) - kurt2  % Ensure kurt2 >= skew2^2 + 1
        ];
        
        % No equality constraints
        ceq = [];
    end  
end



    
    