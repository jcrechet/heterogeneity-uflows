% This script imports in Matlab workspace for quantitative analysis of the model
% empirical statistics computed in Stata
% using CPS-IPUMS data for the U.S. labor force
% OECD.stat cross-country data, and Elsby, Hobijn, Sahin (2013) replication
% files

%% Initialize data structure
data = struct;

%% Age profiles

% Initialize variables.
filename = [data_path, '\CPS_age_profiles.csv'];
opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "age", "e", "u", "ue", "eu", "ee"];
opts.SelectedVariableNames = ["age", "e", "u", "ue", "eu", "ee"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");

% Import the data
data.age = readtable(filename, opts);

% Clear temporary variables
clear opts filename


%% Tenure profiles

% Initialize variables.
filename = [data_path, '\CPS_tenure_profiles.csv'];
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1","tenure", "eu", "ee"];
opts.SelectedVariableNames = ["tenure", "eu", "ee"];
opts.VariableTypes = ["string", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");

% Import the data
data.ten = readtable(filename, opts);

% Clear temporary variables
clear opts filename


%% Unemployment duration profiles

% Initialize variables.
filename = [data_path, '\CPS_udur_profiles.csv'];
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "uduration", "ue"];
opts.SelectedVariableNames = ["uduration", "ue"];
opts.VariableTypes = ["string", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");

% Import the data
data.udur = readtable(filename, opts);

% Clear temporary variables
clear opts filename


%% Compute statistics for calibration

% tmp vectors
u_ = data.age.u;
e_ = data.age.e; 
ue_ = data.age.ue;
eu_ = data.age.eu;
ee_ = data.age.ee;
age_ = data.age.age;

% aggregate rates, age-adjusted, 20-60
age_range = parameters.age_range;
ind_tmp = ( age_ >= age_range(1) & age_ <= age_range(2) );
ue = ( u_(ind_tmp) / sum(u_((ind_tmp))) )'*ue_(ind_tmp);
eu = ( e_(ind_tmp) / sum(e_((ind_tmp))) )'*eu_(ind_tmp);
ee = ( e_(ind_tmp) / sum(e_((ind_tmp))) )'*ee_(ind_tmp);

% ue and eu rate, 20-24
ind_tmp = ( age_ >= 20 & age_ <= 24 );
ue_2024 = ( u_(ind_tmp) / sum(u_((ind_tmp))) )'*ue_(ind_tmp);
eu_2024 = ( e_(ind_tmp) / sum(e_((ind_tmp))) )'*eu_(ind_tmp);

% ue and eu rate, 25-29
ind_tmp = ( age_ >= 25 & age_ <= 29 );
ue_2529 = ( u_(ind_tmp) / sum(u_((ind_tmp))) )'*ue_(ind_tmp);
eu_2529 = ( e_(ind_tmp) / sum(e_((ind_tmp))) )'*eu_(ind_tmp);

% ue and eu rate, 30-39
ind_tmp = ( age_ >= 30 & age_ <= 39 );
ue_3039 = ( u_(ind_tmp) / sum(u_((ind_tmp))) )'*ue_(ind_tmp);
eu_3039 = ( e_(ind_tmp) / sum(e_((ind_tmp))) )'*eu_(ind_tmp);

% ue and eu rate, 40-49
ind_tmp = ( age_ >= 40 & age_ <= 49 );
ue_4049 = ( u_(ind_tmp) / sum(u_((ind_tmp))) )'*ue_(ind_tmp);
eu_4049 = ( e_(ind_tmp) / sum(e_((ind_tmp))) )'*eu_(ind_tmp);

% ue and eu rate, 50-54
ind_tmp = ( age_ >= 50 & age_ <= 54 );
ue_5054 = ( u_(ind_tmp) / sum(u_((ind_tmp))) )'*ue_(ind_tmp);
eu_5054 = ( e_(ind_tmp) / sum(e_((ind_tmp))) )'*eu_(ind_tmp);

% table
calibration_targets = table( ue, eu, ee, ue_2024, ue_2529, ue_3039, ue_4049, ue_5054, eu_2024, eu_2529, eu_3039, eu_4049, eu_5054 );

% log wage variance (see Stata program 2_wage output)
calibration_targets.var_w = .3180192;

% log-wage variance change 1970-1990 (Kambourov Manovskii, 2009)
calibration_targets.var_w_change = 0.57; 

% cleanup
clearvars ue eu ee u_ e_ ue_ eu_ ee_ age_ ue_2024 ue_2529 ue_3039 ue_4049 ue_5054 eu_2024 eu_2529 eu_3039 eu_4049 eu_5054 ind_tmp age_range


%% Cross-country data

filename = [data_path, '\cross_country_data.csv'];

opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["ctry", "country", "ue", "eu", "uemployment_benefits", "ssc_employee", "ssc_employer"];
opts.VariableTypes = ["string", "string", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["ctry", "country"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["ctry", "country"], "EmptyFieldRule", "auto");

% Import the data
data.cross_country = readtable(filename, opts);

clearvars filename opts

%% Cross-country calibration targets

index_US = ( data.cross_country.ctry == "USA" );
index_Eur = ( data.cross_country.ctry == "FRA" | data.cross_country.ctry == "DEU" | data.cross_country.ctry == "ITA" ...
                                      | data.cross_country.ctry == "PRT" | data.cross_country.ctry == "ESP" );

% flows US vs. Europe
calibration_targets.ue_US = data.cross_country.ue( index_US );
calibration_targets.ue_Eur = mean( data.cross_country.ue( index_Eur ) );

calibration_targets.eu_US = data.cross_country.eu( index_US );
calibration_targets.eu_Eur = mean( data.cross_country.eu( index_Eur ) );


% Policy calibration targets

% firing costs
calibration_targets.f_w_Eur = 1;

% unemployment replacement ratio
calibration_targets.b_w_Eur = mean( data.cross_country.uemployment_benefits( index_Eur & data.cross_country.ctry ~= "PRT" ) ) / 100;
calibration_targets.b_w_US  = data.cross_country.uemployment_benefits( index_US ) / 100;

% ss taxes
calibration_targets.tax_US = ( data.cross_country.ssc_employee( index_US ) + data.cross_country.ssc_employer( index_US ) ) / 100;
calibration_targets.tax_Eur = mean( data.cross_country.ssc_employee( index_Eur )  +  data.cross_country.ssc_employer( index_Eur ) ) / 100;

data.calibration_targets = calibration_targets;

clearvars calibration_targets index_Eur index_US





