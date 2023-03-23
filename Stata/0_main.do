* program name: "0_main.do"
* Project: Heterogeneity and unemployment flows across countries
* Jonathan Créchet (2022)
* Created: 2020
* Last updated: March 2023

* Main file for computing CPS statistics used in model calibration/validation

*****
*****

clear
set more off

* path to folder "Replication_files" [fill in for replication]
global mypath = ""

* set directory paths
global working_dir 	    = "$mypath\Stata"
global R_directory		= "$mypath\R"
global Matlab_directory	= "$mypath\Matlab"
global results_path 	= "$mypath\_Results"

******
******

cd "$working_dir"

* 1. Compute worker flows/transition rates
do 1_worker_flows

* 2. Compute wage statistics
do 2_wages

* 3. Cross-country statistics
do 3_cross_country

******
******

