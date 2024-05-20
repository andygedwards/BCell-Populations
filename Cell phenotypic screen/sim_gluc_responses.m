% Code last modified: Roshni (11 August 2023)
% Comments added 2 November 2023

clear all; close all; clc

restoredefaultpath;
addpath('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/Run_Genetic_Algorithm/V_clamp'); % add all subfolders of the working folder's parent folder (this adds all model-specific folders to Matlab's available directories)
addpath('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/Run_Genetic_Algorithm/PoM')
addpath('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/SKNM/Riz_2014_INa_low')
addpath('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/Run_Genetic_Algorithm/Create_optimized_population')

load_par_dir = '/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/Run_Genetic_Algorithm/Create_optimized_population';
nTrials = 3000;
frac_high = 1.00;
exp_name = 'Same';
%%
%Final populations
fname = sprintf('/modParam_scaling_%.0f_%.2f%.mat',nTrials,frac_high);


savedir = sprintf('/Users/Andy/Dropbox/Simula/Research/Current Projects/Pancreatic Islets/Pancreatic-Beta-Cell-Network-SIMULA-main/Cell phenotypic screen/data/%s/%.0f_%.2f%',exp_name,nTrials,frac_high);
if ~ isfolder(savedir)
    mkdir(savedir);
    mkdir(strcat(savedir,'/cells'));
    mkdir(strcat(savedir,'/metadata'));
elseif ~ isfolder(strcat(savedir,'/cells'))
    mkdir(strcat(savedir,'/cells'));
elseif ~ isfolder(strcat(savedir,'/metadata'))
    mkdir(strcat(savedir,'/metadata'));
end

load(sprintf('%s/%s/%.0f_%.2f/%s',load_par_dir,exp_name,nTrials,frac_high,fname))
modParam_names = {'V_GK_max', 'K_GK', 'V_PFK_max', 'K_PFK', 'h_PFK', 'V_GAPDH_max', 'g_Kv', 'g_BK', 'g_CaL', 'g_CaPQ', 'g_CaT', 'g_KATP_hat', 'g_HERG', 'V_hNa', 'n_hNa', 'g_Na', 'V_mNa', 'V_hNa_low','n_hNa_low', 'g_Na_low', 'V_mNa_low'};
modParam_names = convertCharsToStrings(modParam_names);
nModParams = 21;


%% Define Models to Analyze
model = cell(0);
model = [model, {{@Riz2014_init_parameters_INa_low_2, @Riz2014_init_states_INa_low, @Riz2014_rhs_INa_low, 'Beta cell'}}];

model_paramFun = model{1}{1};
model_stateFun = model{1}{2};
model_rhsFun   = model{1}{3};
model_name     = model{1}{4};

%%
[param_vals_baseline, param_names] = model_paramFun();
[param_vals_scaled, modParam_baseline, modParam_vals, modParam_inds] = modifyParams(...
    param_vals_baseline, param_names, modParam_scaling,...
    modParam_names, nTrials);

% Retrieve Baseline State Values
[state_vals_baseline, state_names] = model_stateFun();
%state_vals_baseline has all the inital conditions
nStates = numel(state_vals_baseline);
% EE = zeros(nTrials, 1); TE = zeros(nTrials, 1); 
% peak_INa = zeros(nTrials, 1); v_half_act = zeros(nTrials, 1); vha_r2 = zeros(nTrials, 1); peak_ICa = zeros(nTrials, 1); late_ICa = zeros(nTrials, 1);
% v_half = zeros(nTrials, 1); vh_r2 = zeros(nTrials, 1);

% Compute glucose responses of the population

 if isempty(gcp('nocreate')) % No parallel pool running 
     parpool; % Startup parallel processing pool
 end

 runTime = tic;
 fprintf("Running simulation\n");
 parfor iTrial = 1:nTrials 
    states = state_vals_baseline';
    init_Tstop = 60000;
    param = param_vals_scaled(iTrial,:);
    param(find(strcmp(param_names, 'G'))) = 2;
    options = []; % solver options
    
    % initialization run
    [T_init, Y_init] = ode15s(model_rhsFun, [0, init_Tstop], states, options, param);

    Tstop = 1000000;  % simulation time

    % Low glucose run
    [T_low, Y_low] = ode15s(model_rhsFun, [0, Tstop], Y_init(end,:), options, param);
    fname_low = strcat(savedir,'/cells/',num2str(iTrial),'_low.mat');
    parsave(fname_low,T_low,Y_low,param)

    % High glucose run
    param(find(strcmp(param_names, 'G'))) = 20;
    [T_high, Y_high] = ode15s(model_rhsFun, [0, Tstop], Y_init(end,:), options, param);
    fname_high = strcat(savedir,'/cells/',num2str(iTrial),'_high.mat');
    parsave(fname_high,T_high,Y_high,param)
end

save(strcat(savedir,'/metadata/metadata.mat'),'modParam_names','nModParams','nTrials','model','modParam_inds');


function parsave(fname,T,Y, param)
  save(fname,'T', 'Y','param');
end



