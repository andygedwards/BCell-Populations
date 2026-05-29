% Code last modified: Andy (20 May 2024)

clear; clc; close all

currentDir = pwd;
addpath(genpath(currentDir(1:find(currentDir==filesep,1,'last')-1))); % add all subfolders of the working folder's parent folder (this adds all model-specific folders to Matlab's available directories)

%%
runTime = tic;
fprintf("Running GA simulation\n");
nTrials = 300;

run_ga(nTrials)

toc;

%% run_ga: function description
function [outputs] = run_ga(nTrials)
    input_vector = [0.49, 0.49, 0.15, 0.15, 0.15, 0.15, 0.15];

    %Note gNa stdev(9) replaces later
    initial_population = input_vector;
    nval = 7; % gNa, gNalow, VmNA, VhNa, VmNaLow, VhNaLow and the parameter Frachigh
    ub = [0.49*2, 0.49*2, 0.15*2, 0.15*2, 0.15*2, 0.15*2, 0.50];
    lb = [0, 0, 0, 0, 0, 0, 0.02];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = [];

    options = optimoptions('ga', 'PlotFcn','gaplotbestf', 'UseParallel',true, 'PopulationSize', 25,'MaxGenerations', 20, 'FitnessLimit',100, ...
            'InitialPopulationMatrix', initial_population);
    
    Fitness_Function = @(x)Cost_Function_na(x, nTrials);
    
    [x,fval,exitflag,output, population, costs] = ga(Fitness_Function,nval,A,b,Aeq,beq,lb,ub,nonlcon,options);

filename = sprintf('normal_population_seed_Na.mat');

save (filename, 'population', 'costs', 'x'); %standard_deviations for last generation or stopping criteria of fitness limit reached

outputs = x; %standard_deviations for last generation or stopping criteria when fitness limit reached

end
