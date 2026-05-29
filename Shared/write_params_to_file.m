clear; close all; clc
restoredefaultpath;
home = pwd;
addpath(genpath([home,filesep,'Create Population'])) % Membrane models
addpath(genpath([home,filesep,'Shared'])) % Membrane models
addpath(genpath([home,filesep,'Baseline models'])) % Membrane models

Experiment = 'Combined';
nTrials = 3000;
frac_high = [0.370];
ka = 0.00004;

load_dir = [home,filesep,'Create Population'];

[default_params,default_param_names] = Riz2014_init_parameters_INa_low;

for i = 1:length(frac_high)
    scaling_fname = sprintf('modParam_scaling_%.0f_%.3f.mat',nTrials,frac_high(i)); 
    names_fname = sprintf('modParam_names_%.0f_%.3f.mat',nTrials,frac_high(i));
    load(strcat(load_dir,filesep,Experiment,filesep,sprintf('%.0f',nTrials),filesep,sprintf('%.3f',frac_high(i)),filesep,scaling_fname))
    load(strcat(load_dir,filesep,Experiment,filesep,sprintf('%.0f',nTrials),filesep,sprintf('%.3f',frac_high(i)),filesep,names_fname))
    save_dir = [home,filesep,'Islet Simulator',filesep,'Parameters',filesep,Experiment,filesep,sprintf('%.0f',nTrials),filesep,sprintf('%.3f',frac_high(i))];
    if ~isfolder(save_dir)
        mkdir(save_dir)
    end
    for j = 1:length(modParam_names)
        save_fname = [save_dir,filesep,modParam_names{j},'.txt'];
        save_vector = modParam_scaling(:,j).*(default_params(find(strcmp(default_param_names,modParam_names{j}))));
        writematrix(save_vector, save_fname, 'Delimiter', 'tab'); % Write as a column
    end
end