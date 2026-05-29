% Code last modified: Andy Edwards (October 2024)

clear; clc; close all

currentDir = pwd;
addpath(genpath(currentDir(1:find(currentDir==filesep,1,'last')-1))); % add all subfolders of the working folder's parent folder

%%
runTime = tic;
fprintf("Running PSO simulation\n");
nTrials = 300;

run_pso(nTrials)

toc;

%% run_de: Differential Evolution function
function [outputs] = run_pso(nTrials)
    input_vector = [0.49, 0.49, 0.14, 0.14, 0.14, 0.14, 0.14];

    % Initial values
    nval = 7; % Number of variables: gNa, gNalow, VmNA, VhNa, VmNaLow, VhNaLow, and the parameter Frachigh
    ub = [0.49*2, 0.49*2, 0.15*2, 0.15*2, 0.15*2, 0.15*2, 0.50]; % Upper bounds
    lb = [0, 0, 0, 0, 0, 0, 0.02]; % Lower bounds

    % Setting DE options
    options = optimoptions('particleswarm', 'PlotFcn','pswplotbestf', 'UseParallel',true, ...
                           'SwarmSize', 25, 'MaxIterations', 20, 'FunctionTolerance', 1e-6, ...
                           'InitialSwarmMatrix', input_vector);
    
    % Define the fitness function
    Fitness_Function = @(x)Cost_Function_na(x, nTrials);

    % Run the DE (using particleswarm for DE-like optimization)
    [x,fval,exitflag,output] = particleswarm(Fitness_Function, nval, lb, ub, options);

    % Save results
    filename = sprintf('normal_population_seed_DE_Na.mat');
    save (filename, 'x', 'fval', 'exitflag', 'output'); % Save best solution

    outputs = x; % Best solution found
end
