% Code last modified: Andy Edwards (October 2024)

clear; clc; close all
restoredefaultpath;

home = pwd;
addpath(genpath(home)); % add all subfolders of the global directory

%%
runTime = tic;
fprintf("Running DE simulation\n");

% Set title
title = 'DE optimization 1: Na currents';

% Cost function
ObFunc = @Cost_Function_na;

% Initialize guess for parameter variation, and size of population for fitting
% Order is: [g_Na g_Na_low V_mNa V_hNa V_mNa_low V_hNa_low Frachigh]
initial_vector = [0.49, 0.49, 0.20, 0.20, 0.20, 0.20, 0.50]; % number of decimal places determines the quantization of the search function (10-fold)
nTrials = 200;

% run the DE
[best, OptErr, OptParams, OptIters, OptFname] = run_DE (nTrials, initial_vector, ObFunc, title);

toc;

function [bestmem, bestval, bestFctParams, nrOfIterations, resultFileName] = run_DE(n, init, func, title)

    % Set title
    optimInfo.title = title;
    
    % Specify objective function
    objFctHandle = func;

    for i = 1: length(init)
        temp = regexp(num2str(init(i)),'\.','split');
        len(i) = length(temp{2});
        quant = 10^-(max(len)+1);
    end

    % 1. column: parameter names
    % 2. column: parameter ranges
    % 3. column: parameter quantizations
    paramDefCell = {
	    'g_Na', [init(1)/2 init(1)*2], quant;
        'g_Na_low', [init(2)/2 init(2)*2], quant; 
	    'V_mNa', [init(3)/3 init(3)*3], quant; 
        'V_hNa', [init(4)/3 init(4)*3], quant; 
        'V_mNa_low', [init(5)/3 init(5)*3], quant; 
        'V_hNa_low', [init(6)/3 init(6)*3], quant; 
        'Frachigh', [init(7)/5 init(7)*2], quant; 
    };
    
    % Set initial parameter values in struct objFctParams 
    objFctParams.g_Na = init(1);
    objFctParams.g_Na_low = init(2);
    objFctParams.V_mNa = init(3);
    objFctParams.V_hNa = init(4);
    objFctParams.V_mNa_low = init(5);
    objFctParams.V_hNa_low = init(6);
    objFctParams.Frachigh = init(7);
    
    objFctSettings = n; %(nTrials)

    DEParams = getdefaultparams;
    DEParams.NP = 75*length(init);

    % Not parallelized
    DEParams.feedSlaveProc = 0;

    % Set times
    DEParams.maxiter = 1;
    DEParams.maxtime = inf; % in seconds
    DEParams.maxclock = [];
    DEParams.minimizeValue = 1;

    % Display options
    DEParams.infoIterations = 1;
    DEParams.infoPeriod = 60;  % in seconds

    % Do not send E-mails
    emailParams = [];

    % Set random state in order to always use the same population members here
    setrandomseed(0);
    % Start differential evolution
    [bestmem, bestval, bestFctParams, nrOfIterations, resultFileName] = differentialevolution(...
	DEParams, paramDefCell, objFctHandle, objFctSettings, objFctParams, emailParams, optimInfo);

    disp(' ');
    disp('Best parameter set returned by function differential evolution:');
    disp(bestFctParams);
end




