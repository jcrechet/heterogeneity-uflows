% Project: Heterogeneity in labor mobility and unemployment flows across countries
% Créchet, 2023
% created 2020 - last updated March 2023

% Main Matlab script for calibration and analysis of quantitative model

%% 1. Preliminaries

clear, clc
timer = tic;

% set path to Replication folders
mypath = '';

% working directory path
dir_path = [mypath,'\Matlab\Full_model'];
% data path
data_path = [mypath,'\Matlab\Data'];
% results path
results_path = [mypath,'\_Results'];

% Set current folder
cd(dir_path)

% Add functions to working path
addpath('Functions');

% Call script that creates parameter structure
run('SetParameters.m')

% Call script to import empirical statistics to Matlab Workspace
run('ImportData.m')


%% 2. QUANTITATIVE ANALYSIS

%% 2.1. Benchmark model outcomes

% evaluate model at benchmark param values
parameters.display_equil_iterations = 1;
parameters.equilibrium = 1;
parameters.calibration = 0;
[~, moments, variables, ~] = equilibrium([], parameters, []);

% display model select outcomes
disp('************')
disp('model select outcomes:')
disp(moments.agg)

% save outcomes of benchmark calibrated model
parameters.display_equil_iterations = 0;
save("Results/Benchmark.mat", "variables", "moments", "parameters", '-append')

% calibration outcomes
calibration_outcomes(data, parameters, moments, results_path)


%% 2.2. Elasticity decompositions
run("Decompositions.m")


%% 2.3. Calibration of alternative models
run("Alternative_models.m")


%% 2.4. Accounting for cross-country differences across alternative models
run("Accounting.m")


%% 2.5. Secular US-Europe u-differences
run('Secular_urate.m')

%%
disp('total time:')
total_time = toc(timer);
disp(total_time);
