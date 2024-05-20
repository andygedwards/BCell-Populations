function Total_Error = Cost_Function_exo_cal(input_vec, nTrials)

% Code last modified: Roshni (11 August 2023)
% Comments added 2 November 2023

% Fraction of cells expressing high voltage INa
frac_high = 0.15;   %15 percent INa high
frac_low = 1-frac_high; %85% INa low

st_dev_vector = zeros(1, 22);
st_dev_vector(10) = input_vec(1);  %gCaL
st_dev_vector(11) = input_vec(2);  %gCaPQ
st_dev_vector(19) = input_vec(3);  %V_hNa_low

modParam_stdev_noNa = st_dev_vector(1:14); % standard deviations (normalized to the mean) of the parameter names below
INa_high_modParam_stdev = st_dev_vector(15:18);
INa_low_modParam_stdev = st_dev_vector(19:22);


% modParam_stdev_noNa = [0.71 0.90 0.03 0.09 0.14 0.27 0.37 0.80 0.49 0.34 0.40 0.42 0.89 0.30]; % standard deviations (normalized to the mean) of the parameter names below
% INa_high_modParam_stdev = [0.15, 0.15, 0.49, 0.15];
% INa_low_modParam_stdev = [0.15, 0.15, 0.49, 0.15];

%Note gNa stdev(9) replaced later in code
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
INa_high_modParam_scaling = getScalingFactors_INa_inact(INa_high_modParam_stdev, INa_high_modParam_names, frac_high*nTrials);
INa_high_modParam_scaling = [INa_high_modParam_scaling;zeros(frac_low*nTrials,INa_high_nModParams)];

INa_low_modParam_names = {'V_hNa_low','n_hNa_low', 'g_Na_low', 'V_mNa_low'};
INa_low_nModParams = numel(INa_low_modParam_names);
INa_low_modParam_scaling = getScalingFactors_INa_inact(INa_low_modParam_stdev, INa_low_modParam_names, frac_low*nTrials);
INa_low_modParam_scaling = [zeros(frac_high*nTrials,INa_low_nModParams);INa_low_modParam_scaling];

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
%state_vals_baseline has al l the inital conditions for 15 isletparams
nStates = numel(state_vals_baseline);
EE = zeros(nTrials, 1); TE = zeros(nTrials, 1); 
peak_INa = zeros(nTrials, 1); v_half_act = zeros(nTrials, 1); vha_r2 = zeros(nTrials, 1); peak_ICa = zeros(nTrials, 1); late_ICa = zeros(nTrials, 1);
v_half = zeros(nTrials, 1); vh_r2 = zeros(nTrials, 1);

% Compute the output metrics for the population
try
for iTrial = 1:nTrials 
    state_vals_output = state_vals_baseline';
    
    %[EE, LE, TE, v_half, vh_r2, peak_INa, v_half_act, vha_r2, peak_ICa, late_ICa, ramp_IK, Peak_IK_minus20, Peak_IK_max, v_half_act_IK, n_act_IK] = run_v_clamp(param_vals_scaled(iTrial,:), param_names,...
         %state_vals_output, state_names, model_rhsFun);
    [EE(iTrial, 1), ~, TE(iTrial,1)] = Exocytosis(param_vals_scaled(iTrial,:), param_names, state_vals_output, state_names, model_rhsFun);
    [peak_INa(iTrial, 1), v_half_act(iTrial, 1), vha_r2(iTrial, 1), peak_ICa(iTrial, 1), late_ICa(iTrial, 1)] = PeakCurrents(param_vals_scaled(iTrial,:), param_names, state_vals_output, state_names, model_rhsFun);
    [v_half(iTrial, 1), vh_r2(iTrial, 1)] = INa_inact(param_vals_scaled(iTrial,:), param_names, state_vals_output, state_names, model_rhsFun);

end

%%  EE
%Simulated distribution characteristics
mu_EE_trial = mean(EE); sigma_EE_trial = std(EE);
skewness_EE_trial = skewness(EE); kurtosis_EE_trial = kurtosis(EE);
 
%Experimental distribution characteristics from Camunas-Soler (2019)
mu_EE_dist = 3.88; sigma_EE_dist = 6.56;
skewness_EE_dist = 4.78; kurtosis_EE_dist = 31.96;

EE_error1 = abs((sigma_EE_trial/mu_EE_trial - sigma_EE_dist/mu_EE_dist)/ (sigma_EE_dist/mu_EE_dist) );
EE_error2 = abs((skewness_EE_trial - skewness_EE_dist)/skewness_EE_dist);
EE_error3 = abs((kurtosis_EE_trial - kurtosis_EE_dist)/kurtosis_EE_dist);

EE_Total_Error = 0.5*100^2 * (EE_error1 + EE_error2 + EE_error3);

%% TE
%Simulated distribution characteristics
mu_TE_trial = mean(TE); sigma_TE_trial = std(TE);
skewness_TE_trial = skewness(TE); kurtosis_TE_trial = kurtosis(TE);

%Experimental distribution characteristics from Camunas-Soler (2019)
mu_TE_dist = 16.78; sigma_TE_dist = 24.43;
skewness_TE_dist = 5.11; kurtosis_TE_dist = 40.05;

TE_error1 = abs((sigma_TE_trial/mu_TE_trial - sigma_TE_dist/mu_TE_dist)/ (sigma_TE_dist/mu_TE_dist) );
TE_error2 = abs((skewness_TE_trial - skewness_TE_dist)/skewness_TE_dist);
TE_error3 = abs((kurtosis_TE_trial - kurtosis_TE_dist)/kurtosis_TE_dist);

TE_Total_Error = 0.5*100^2 * (TE_error1 + TE_error2 + TE_error3);


%% peak ICa
%Simulated distribution characteristics
mu_peak_ICa_trial = mean(peak_ICa); sigma_peak_ICa_trial = std(peak_ICa);
skewness_peak_ICa_trial = skewness(peak_ICa); kurtosis_peak_ICa_trial = kurtosis(peak_ICa);

%Experimental distribution characteristics from Camunas-Soler (2019)
mu_peak_ICa_dist = -3.21; sigma_peak_ICa_dist = 2.46;
skewness_peak_ICa_dist =  -1.24; kurtosis_peak_ICa_dist = 5.09;

peak_ICa_error1 = abs((sigma_peak_ICa_trial/mu_peak_ICa_trial - sigma_peak_ICa_dist/mu_peak_ICa_dist)/ (sigma_peak_ICa_dist/mu_peak_ICa_dist) );
peak_ICa_error2 = abs((skewness_peak_ICa_trial - skewness_peak_ICa_dist)/skewness_peak_ICa_dist);
peak_ICa_error3 = abs((kurtosis_peak_ICa_trial - kurtosis_peak_ICa_dist)/kurtosis_peak_ICa_dist);

peak_ICa_Total_Error = 0.5*100^2 * (peak_ICa_error1 + peak_ICa_error2 + peak_ICa_error3);

%% Late ICaL
%Simulated distribution characteristics
mu_late_ICa_trial = mean(late_ICa); sigma_late_ICa_trial = std(late_ICa);
skewness_late_ICa_trial = skewness(late_ICa); kurtosis_late_ICa_trial = kurtosis(late_ICa);

%Experimental distribution characteristics from Camunas-Soler (2019)
mu_late_ICa_dist = -1.04; sigma_late_ICa_dist = 0.99;
skewness_late_ICa_dist =  -1.14; kurtosis_late_ICa_dist = 4.03;

late_ICa_error1 = abs((sigma_late_ICa_trial/mu_late_ICa_trial - sigma_late_ICa_dist/mu_late_ICa_dist)/ (sigma_late_ICa_dist/mu_late_ICa_dist) );
late_ICa_error2 = abs((skewness_late_ICa_trial - skewness_late_ICa_dist)/skewness_late_ICa_dist);
late_ICa_error3 = abs((kurtosis_late_ICa_trial - kurtosis_late_ICa_dist)/kurtosis_late_ICa_dist);

late_ICa_Total_Error = 0.5*100^2 * (late_ICa_error1 + late_ICa_error2 + late_ICa_error3);


%%
% GA minimizes the Total_Error
Total_Error = EE_Total_Error + TE_Total_Error + peak_ICa_Total_Error + late_ICa_Total_Error;

if Total_Error < 18000
    parsave2(Total_Error, modParam_scaling, modParam_stdev, peak_INa, v_half_act, peak_ICa, late_ICa, ...
        EE, TE, v_half,  ....
        mu_EE_trial, sigma_EE_trial, skewness_EE_trial, kurtosis_EE_trial, ...
        mu_TE_trial, sigma_TE_trial, skewness_TE_trial, kurtosis_TE_trial,...
        mu_late_ICa_trial, sigma_late_ICa_trial, skewness_late_ICa_trial, kurtosis_late_ICa_trial,...
        mu_peak_ICa_trial, sigma_peak_ICa_trial, skewness_peak_ICa_trial, kurtosis_peak_ICa_trial)

end


catch
    disp('Error encountered');
    Total_Error = 10e7;
end


function parsave2(Total_Error, modParam_scaling, modParam_stdev, peak_INa, v_half_act, peak_ICa, late_ICa, ...
        EE, TE, v_half,  ....
        mu_EE_trial, sigma_EE_trial, skewness_EE_trial, kurtosis_EE_trial, ...
        mu_TE_trial, sigma_TE_trial, skewness_TE_trial, kurtosis_TE_trial,...
        mu_late_ICa_trial, sigma_late_ICa_trial, skewness_late_ICa_trial, kurtosis_late_ICa_trial,...
        mu_peak_ICa_trial, sigma_peak_ICa_trial, skewness_peak_ICa_trial, kurtosis_peak_ICa_trial)

h = figure(); histogram(EE, 25, 'Normalization', 'probability')
heading = sprintf('Error%.2f mu%.2f sigma%.2f skew%.2f kur%.2f', Total_Error, mu_EE_trial, sigma_EE_trial, skewness_EE_trial, kurtosis_EE_trial);
title(heading)
saveas(h,sprintf('FIG_EE%.2f.png',Total_Error)); 

h1 = figure(); histogram(TE, 25, 'Normalization', 'probability')
heading1 = sprintf('Error%.2f mu%.2f sigma%.2f skew%.2f kur%.2f', Total_Error, mu_TE_trial, sigma_TE_trial, skewness_TE_trial, kurtosis_TE_trial);
title(heading1)
saveas(h1,sprintf('FIG_TE%.2f.png',Total_Error)); 

h2 = figure(); histogram(peak_ICa, 25, 'Normalization', 'probability')
heading2 = sprintf('Error%.2f mu%.2f sigma%.2f skew%.2f kur%.2f', Total_Error, mu_peak_ICa_trial, sigma_peak_ICa_trial, skewness_peak_ICa_trial, kurtosis_peak_ICa_trial);
title(heading2)
saveas(h2,sprintf('FIG_earlyCa%.2f.png',Total_Error)); 

h3 = figure(); histogram(late_ICa, 25, 'Normalization', 'probability')
heading3 = sprintf('Error%.2f mu%.2f sigma%.2f skew%.2f kur%.2f', Total_Error, mu_late_ICa_trial, sigma_late_ICa_trial, skewness_late_ICa_trial, kurtosis_late_ICa_trial);
title(heading3)
saveas(h3,sprintf('FIG_lateCa%.2f.png',Total_Error)); 


scaling_name = sprintf('modParam_scaling_%.2f.mat', Total_Error); save(scaling_name, 'modParam_scaling');
stdev_name = sprintf('modParam_stdev_%.2f.mat', Total_Error); save(stdev_name, 'modParam_stdev');
peak_INa_name = sprintf('peak_INa_%.2f.mat', Total_Error); save(peak_INa_name, 'peak_INa');
v_half_act_name = sprintf('v_half_act_%.2f.mat', Total_Error); save(v_half_act_name, 'v_half_act');
peak_ICa_name = sprintf('peak_ICa_%.2f.mat', Total_Error); save(peak_ICa_name, 'peak_ICa');
late_ICa_name = sprintf('late_ICa_%.2f.mat', Total_Error); save(late_ICa_name, 'late_ICa');
EE_name = sprintf('EE_%.2f.mat', Total_Error); save(EE_name, 'EE');
TE_name = sprintf('TE_%.2f.mat', Total_Error); save(TE_name, 'TE');
v_half_name = sprintf('v_half_%.2f.mat', Total_Error); save(v_half_name, 'v_half');

end
  
end



    
    