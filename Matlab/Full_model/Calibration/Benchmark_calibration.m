% Project: Heterogeneity in labor mobility and unemployment flows across countries
% Créchet, 2022
% Program name: "Main.m"
% created 2020 - updated September 2022

% Script for calibration of benchmark model

%% 1. preliminaries

clear, clc

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


%% 3. calibration algorithm

%% step 1: random grid search

% seed state
rng(1)

% options
parameters.calibration = true;
run_rand_search = true;
run_ga_minimization = false;

% run function
calibration_benchmark(parameters, parameters.x0, run_rand_search, run_ga_minimization, []);


%% step 2: ga minization (combined with Matlab function pattersearch)

run_rand_search = false;
run_ga_minimization = true;
run_ga_initial = true;

calibration_benchmark(parameters, [], run_rand_search, run_ga_minimization, run_ga_initial);

% note: step 1 must be run at least one before step 2. Setting
% run_ga_initial = false reevaluate the objective function over the last best grid of parameters in memory. Use this option 
% if you want to start from a preexisting grid. Running step 1 overwrite
% the preexisting grid. Run step 1 to start from new parameter points.
