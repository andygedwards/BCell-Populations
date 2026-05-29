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
    input_vector = [0.34, 0.4]; % gCaPQ, gCaL

    %Note gNa stdev(9) replaced later in code
    initial_population = input_vector;
    nval = 3; % gCaPQ, gCaL, V_hNa_low
    ub = input_vector.*2;
    lb = zeros(3, 1);
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = [];

    options = optimoptions('ga', 'PlotFcn','gaplotbestf', 'UseParallel',true, 'PopulationSize', 25,'MaxGenerations', 20, 'FitnessLimit',500, ...
            'InitialPopulationMatrix', initial_population);
    
    Fitness_Function = @(x)Cost_Function_exo_cal(x, nTrials);
    
    [x,fval,exitflag,output, population, costs] = ga(Fitness_Function,nval,A,b,Aeq,beq,lb,ub,nonlcon,options);

filename = sprintf('normal_population_seed_exo_caL.mat');

save (filename, 'population', 'costs', 'x');

outputs = x; %standard_deviations for last generation or stopping criteria when fitness limit reached
end
