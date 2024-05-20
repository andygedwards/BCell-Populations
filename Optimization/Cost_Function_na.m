function Total_Error = Cost_Function_na(input_vec, nTrials)

    % Code last modified: Roshni (11 August 2023)
    % Comments added 2 November 2023
    
    % Fraction of cells expressing high voltage INa
    frac_high = round(input_vec(7),2);  %15 percent INa high
    frac_low = round(1-frac_high,2);   %85% INa low
    
    trials_low = round(nTrials*frac_low);
    trials_high = nTrials - trials_low;
    
    st_dev_vector = zeros(1, 22);
    
    st_dev_vector(17) = input_vec(1); %gNa
    st_dev_vector(21) = input_vec(2); %gNa_low
    st_dev_vector(15) = input_vec(3); %V_hNa
    st_dev_vector(19) = input_vec(4); %V_hNa_low
    st_dev_vector(18) = input_vec(5); %V_mN
    st_dev_vector(22) = input_vec(6); %V_mN_low
    
    modParam_stdev_noNa = st_dev_vector(1:14); % standard deviations (normalized to the mean) of the parameter names below
    INa_high_modParam_stdev = st_dev_vector(15:18);
    INa_low_modParam_stdev = st_dev_vector(19:22);
    
    % modParam_stdev_noNa = [0.71 0.90 0.03 0.09 0.14 0.27 0.37 0.80 0.49 0.34 0.40 0.42 0.89 0.30]; % standard deviations (normalized to the mean) of the parameter names below
    % INa_high_modParam_stdev = [0.15, 0.15, 0.49, 0.15];
    % INa_low_modParam_stdev = [0.15, 0.15, 0.49, 0.15];
    
    modParam_names_noNa = {'V_GK_max', 'K_GK', 'V_PFK_max', 'K_PFK', 'h_PFK', 'V_GAPDH_max', 'g_Kv', 'g_BK', 'g_Na', 'g_CaL', 'g_CaPQ', 'g_CaT', 'g_KATP_hat', 'g_HERG'};
    nModParams_noNa = numel(modParam_names_noNa);
    modParam_scaling_noNa = getScalingFactors(modParam_stdev_noNa, nModParams_noNa, nTrials);
    
    %% Define Models to Analyze
    model = cell(0);
    model = [model, {{@Riz2014_init_parameters_INa_low, @Riz2014_init_states_INa_low, @Riz2014_rhs_INa_low, 'Beta cell'}}];
    
    
    %% INa_high_low
    modParam_scaling_noNa(:,find(strcmp(modParam_names_noNa,'g_Na'))) = [];
    modParam_stdev_noNa(:,find(strcmp(modParam_names_noNa,'g_Na'))) = [];
    modParam_names_noNa(find(strcmp(modParam_names_noNa,'g_Na'))) = [];
    nModParams = nModParams_noNa-1;
    
    INa_high_modParam_names = {'V_hNa', 'n_hNa', 'g_Na', 'V_mNa'};
    INa_high_nModParams = numel(INa_high_modParam_names);
    INa_high_modParam_scaling = getScalingFactors_INa_inact(INa_high_modParam_stdev, INa_high_modParam_names, trials_high);
    INa_high_modParam_scaling = [INa_high_modParam_scaling; zeros(trials_low,INa_high_nModParams)];
    
    % INa_high_modParam_scaling = getScalingFactors_INa_inact(INa_high_modParam_stdev, INa_high_modParam_names, frac_high*nTrials);
    % INa_high_modParam_scaling = [INa_high_modParam_scaling;zeros(frac_low*nTrials,INa_high_nModParams)];
    
    INa_low_modParam_names = {'V_hNa_low','n_hNa_low', 'g_Na_low', 'V_mNa_low'};
    INa_low_nModParams = numel(INa_low_modParam_names);
    INa_low_modParam_scaling = getScalingFactors_INa_inact(INa_low_modParam_stdev, INa_low_modParam_names, trials_low);
    INa_low_modParam_scaling = [zeros(trials_high,INa_low_nModParams);INa_low_modParam_scaling];
    
    % INa_low_modParam_scaling = getScalingFactors_INa_inact(INa_low_modParam_stdev, INa_low_modParam_names, frac_low*nTrials);
    % INa_low_modParam_scaling = [zeros(frac_high*nTrials,INa_low_nModParams);INa_low_modParam_scaling];
    
    modParam_scaling = [modParam_scaling_noNa,INa_high_modParam_scaling,INa_low_modParam_scaling];
    modParam_names = [modParam_names_noNa,INa_high_modParam_names, INa_low_modParam_names];
    modParam_stdev = [modParam_stdev_noNa,INa_high_modParam_stdev,INa_low_modParam_stdev];
    nModParams = nModParams+INa_high_nModParams+INa_low_nModParams;
    
    % figureFolder_general = [output_general, filesep, 'Figures'];
    % if exist(figureFolder_general, 'dir') ~= 7
    %     mkdir(figureFolder_general);
    % end
    % plotScalingFactors(modParam_scaling, modParam_names, nModParams, figureFolder_general);
    
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
    %state_vals_baseline has all the inital conditions for 15 isletparams
    nStates = numel(state_vals_baseline);
    EE = zeros(nTrials, 1); TE = zeros(nTrials, 1); 
    peak_INa = zeros(nTrials, 1); v_half_act = zeros(nTrials, 1); vha_r2 = zeros(nTrials, 1); peak_ICa = zeros(nTrials, 1); late_ICa = zeros(nTrials, 1);
    v_half = zeros(nTrials, 1); vh_r2 = zeros(nTrials, 1);
    
    % Compute the output metrics for the population
    try
        for iTrial = 10:1000:1010%1:nTrials 
            state_vals_output = state_vals_baseline';
            %[EE, LE, TE, v_half, vh_r2, peak_INa, v_half_act, vha_r2, peak_ICa, late_ICa, ramp_IK, Peak_IK_minus20, Peak_IK_max, v_half_act_IK, n_act_IK] = run_v_clamp(param_vals_scaled(iTrial,:), param_names,...
                 %state_vals_output, state_names, model_rhsFun);
            [EE(iTrial, 1), ~, TE(iTrial,1)] = Exocytosis(param_vals_scaled(iTrial,:), param_names, state_vals_output, state_names, model_rhsFun);
            [peak_INa(iTrial, 1), v_half_act(iTrial, 1), vha_r2(iTrial, 1), peak_ICa(iTrial, 1), late_ICa(iTrial, 1)] = PeakCurrents(param_vals_scaled(iTrial,:), param_names, state_vals_output, state_names, model_rhsFun);
            [v_half(iTrial, 1), vh_r2(iTrial, 1)] = INa_inact(param_vals_scaled(iTrial,:), param_names, state_vals_output, state_names, model_rhsFun);
        end
    
        %% peak INa
        %Simulated distribution characteristics
        mu_peak_INa_trial = mean(peak_INa); sigma_peak_INa_trial = std(peak_INa);
        skewness_peak_INa_trial = skewness(peak_INa); kurtosis_peak_INa_trial = kurtosis(peak_INa);
        
        %Experimental distribution characteristics from Camunas-Soler (2019)
        mu_peak_INa_dist = -13.56; sigma_peak_INa_dist = 15.62;
        skewness_peak_INa_dist = -1.94; kurtosis_peak_INa_dist = 6.94;
        
        peak_INa_error1 = abs((sigma_peak_INa_trial/mu_peak_INa_trial - sigma_peak_INa_dist/mu_peak_INa_dist)/ (sigma_peak_INa_dist/mu_peak_INa_dist) );
        peak_INa_error2 = abs((skewness_peak_INa_trial - skewness_peak_INa_dist)/skewness_peak_INa_dist);
        peak_INa_error3 = abs((kurtosis_peak_INa_trial - kurtosis_peak_INa_dist)/kurtosis_peak_INa_dist);
        
        peak_INa_Total_Error = 0.5*100^2 * (peak_INa_error1 + peak_INa_error2 + peak_INa_error3);
        
        %% v_half 
        %Fit double Gaussian to population
        y = histogram(v_half, 'BinWidth', 2);
        counts = y.BinCounts;   
        edges = y.BinEdges; middle_bins = (edges(2:end) + edges(1:end-1))/2;
        close all
        % Define the x and y data
        xData = middle_bins'; yData = counts';
        doubleGaussianModel = @(a1, mu1, sigma1, a2, mu2, sigma2, x) ...
            a1 * exp(-(x - mu1).^2 / (2 * sigma1^2)) + ...
            a2 * exp(-(x - mu2).^2 / (2 * sigma2^2));
        
        % Initial estimates of a, mu and sigma
        initialGuesses = [1, -45, 4.135, 1, -65, 3.673];
        
        %Simulated distribution characteristics
        fitResult = fit(xData, yData, doubleGaussianModel, 'Start', initialGuesses);
        mu1_fit = fitResult.mu1; mu2_fit = fitResult.mu2;
        sigma1_fit = abs(fitResult.sigma1); sigma2_fit = abs(fitResult.sigma2);
        a1_fit = fitResult.a1; a2_fit = fitResult.a2;
        a1a2_fit = a2_fit/a1_fit;
        
        %Experimental distribution characteristics from Camunas-Soler (2019)
        mu1_expt = -47.15; mu2_expt = -63.7;
        sigma1_expt = 4.74; sigma2_expt = 3.508;
        a1_expt = 9.444; a2_expt = 15.01;
        a1a2_expt = a2_expt/a1_expt;
        
        vhalf_error1= abs((mu1_fit - mu1_expt)/mu1_expt);
        vhalf_error2= abs((mu2_fit - mu2_expt)/mu2_expt);
        vhalf_error3= abs((sigma1_fit - sigma1_expt)/sigma1_expt);
        vhalf_error4= abs((sigma2_fit - sigma2_expt)/sigma2_expt);
        vhalf_error5 = abs((a1a2_fit - a1a2_expt)/a1a2_expt);
        
        vhalf_Total_Error = vhalf_error1 + vhalf_error2 + vhalf_error3 + vhalf_error4 + vhalf_error5;
       
        Total_Error = peak_INa_Total_Error + 2*vhalf_Total_Error;
        
%         if Total_Error < 8000
%             parsave2(Total_Error, modParam_scaling, modParam_stdev, peak_INa, v_half_act, peak_ICa, late_ICa, ...
%                 EE, TE, v_half,  frac_high,....
%                 mu_peak_INa_trial, sigma_peak_INa_trial, skewness_peak_INa_trial, kurtosis_peak_INa_trial, ...
%                 mu1_fit, mu2_fit, sigma1_fit, sigma2_fit)         
%         end
    catch
        disp('Error encountered');
        Total_Error = 10e7;
    end
    
    function parsave2(Total_Error, modParam_scaling, modParam_stdev, peak_INa, v_half_act, peak_ICa, late_ICa, ...
            EE, TE, v_half, frac_high, ....
            mu_peak_INa_trial, sigma_peak_INa_trial, skewness_peak_INa_trial, kurtosis_peak_INa_trial, ...
            mu1_fit, mu2_fit, sigma1_fit, sigma2_fit)
    
        h = figure(); histogram(peak_INa, 25, 'Normalization', 'probability')
        heading = sprintf('Error%.2f mu%.2f sigma%.2f skew%.2f kur%.2f', Total_Error, mu_peak_INa_trial, sigma_peak_INa_trial, skewness_peak_INa_trial, kurtosis_peak_INa_trial);
        title(heading)
        saveas(h,sprintf('FIG_peakINa%.2f.png',Total_Error)); 
        
        h1 = figure(); histogram(v_half, 25, 'Normalization', 'probability')
        heading1 = sprintf('Error%.2f mu1%.2f sigma1%.2f mu2%.2f sigma2%.2f', Total_Error, mu1_fit, mu2_fit, sigma1_fit, sigma2_fit);
        title(heading1)
        saveas(h1,sprintf('FIG_vhalf%.2f.png',Total_Error)); 
        
        scaling_name = sprintf('modParam_scaling_%.2f.mat', Total_Error); save(scaling_name, 'modParam_scaling');
        stdev_name = sprintf('modParam_stdev_%.2f.mat', Total_Error); save(stdev_name, 'modParam_stdev');
        peak_INa_name = sprintf('peak_INa_%.2f.mat', Total_Error); save(peak_INa_name, 'peak_INa');
        v_half_act_name = sprintf('v_half_act_%.2f.mat', Total_Error); save(v_half_act_name, 'v_half_act');
        peak_ICa_name = sprintf('peak_ICa_%.2f.mat', Total_Error); save(peak_ICa_name, 'peak_ICa');
        late_ICa_name = sprintf('late_ICa_%.2f.mat', Total_Error); save(late_ICa_name, 'late_ICa');
        EE_name = sprintf('EE_%.2f.mat', Total_Error); save(EE_name, 'EE');
        TE_name = sprintf('TE_%.2f.mat', Total_Error); save(TE_name, 'TE');
        v_half_name = sprintf('v_half_%.2f.mat', Total_Error); save(v_half_name, 'v_half');
        frac_high_name = sprintf('frac_high_%.2f.mat', Total_Error); save(frac_high_name, 'frac_high');
    end
end



    
    