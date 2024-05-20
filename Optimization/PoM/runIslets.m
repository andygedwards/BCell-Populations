close all; clear; clc; 
currentDir = pwd;
addpath(genpath(currentDir(1:find(currentDir==filesep,1,'last')-1))); % add all subfolders of the working folder's parent folder (this adds all model-specific folders to Matlab's available directories)

%% Simulation output and simulation parameters

%saves the output of the ode solver in this uniquely named folder
addpath(genpath(currentDir(1:find(currentDir==filesep,1,'last')-1)));
output_general = [pwd, filesep, 'output',datestr(now,'yyyymmddTHHMMSS')];
if exist(output_general, 'dir') ~= 7
    mkdir(output_general);
end

%% Quick-Access Analysis Settings
nTrials = 50;          % number of different parameter combinations (target = 3000)
modParam_stdev = [0.71 0.90 0.03 0.09 0.14 0.27 0.37 0.80 0.49 0.34 0.40 0.42 0.89 0.30]; % standard deviations (normalized to the mean) of the parameter names below
modParam_names = {'V_GK_max', 'K_GK', 'V_PFK_max', 'K_PFK', 'h_PFK', 'V_GAPDH_max', 'g_Kv', 'g_BK', 'g_Na', 'g_CaL', 'g_CaPQ', 'g_CaT', 'g_KATP_hat', 'g_HERG'};
              
%% Define Models to Analyze
% model = [model, {{parameter_function_handle, state_function_handle, rhs_function_handle, 'modelName'}}]; % example format (comment-out when running code)

% Individual models (code lines) can be commented-out without affecting how the code runs/works
model = cell(0);
model = [model, {{@Riz2014_init_parameters, @Riz2014_init_states, @Riz2014_rhs, 'islet'}}];
nModels = numel(model);

%% Setup Model Parameters to be Perturbed

cd([pwd, filesep,'POM_Scaling_Factors']);
scalingFactorFile = 'SAVED_scalingFactors.mat';
if isfile(scalingFactorFile)
    load(scalingFactorFile);
    fprintf("Scaling factors loaded.\n");
else
    nModParams = numel(modParam_names);
    modParam_scaling = getScalingFactors(modParam_stdev, nModParams, nTrials);
%     save(scalingFactorFile, 'modParam_scaling', 'modParam_names', 'nModParams', 'modParam_stdev', 'nTrials', 'nStims', 'stim_freq', '-v7.3');
    save(scalingFactorFile, 'modParam_scaling', 'modParam_names', 'nModParams', 'modParam_stdev', 'nTrials', '-v7.3');
    fprintf("Scaling factors generated.\n");
end
cd ..

figureFolder_general = [pwd, filesep, 'Figures'];
if exist(figureFolder_general, 'dir') ~= 7
    mkdir(figureFolder_general);
end
plotScalingFactors(modParam_scaling, modParam_names, nModParams, figureFolder_general);

modParam_vals = cell(nModels, 1);
IC_vals = cell(nModels, 1);
IC_names = cell(nModels, 1);
trial_flags = cell(nModels, 1);

param_vals_scaled = cell(nModels,1); 
modParam_baseline = cell(nModels,1); 

model_paramFun = model{1}{1};
model_stateFun = model{1}{2};
model_rhsFun   = model{1}{3};
model_name     = model{1}{4};

[param_vals_baseline, param_names] = model_paramFun();
[param_vals_scaled, modParam_baseline, modParam_vals] = modifyParams(...
    param_vals_baseline, param_names, modParam_scaling,...
    modParam_names, nTrials);

 % Retrieve Baseline State Values
 [state_vals_baseline, state_names] = model_stateFun();
 %state_vals_baseline has all the inital conditions for 14 islet params
 nStates = numel(state_vals_baseline);
 
 % Simulation Parameters
 Tstop = 1000*10*60;  % simulation time
 options = [];
 time_range = [0, Tstop];
 
 if isempty(gcp('nocreate')) % No parallel pool running 
     parpool; % Startup parallel processing pool
 end
   
 % Solve the model for each parameter trial and store ICs
 runTime = tic;
 fprintf("Running simulation\n");
 parfor iTrial = 1:nTrials  
     state_vals_output = state_vals_baseline';
     [time_output, state_vals_output] = ode15s(model_rhsFun, ...
         time_range, state_vals_output, ...
         options, param_vals_scaled(iTrial,:));
      parsave(output_general,time_output,state_vals_output, iTrial)
 end

fprintf("\n%.3fs\n", toc(runTime));
clear all;

function parsave(output_general,T,Y, i)
  saveoutputfilepath = [output_general, filesep,num2str(i)];
  save(saveoutputfilepath,'T', 'Y');
end






