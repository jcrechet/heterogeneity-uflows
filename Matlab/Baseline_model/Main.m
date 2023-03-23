% Créchet J.
% created: december 2020; last updated March 2023.
% project: Heterogeneity in labor mobility and unemployment across
% countries
% compute semi-elasticities of unemployment with respect to firing costs
% in a DMP model with endogenous separations and heterogeneous
% mobility (stochastic matching)

clear 
clc

tic

% path to Replication folder
mypath = '';

% Working directory path
dir_path = [mypath,'\Matlab\Baseline_model'];

% Result paths
results_path = [mypath,'\_Results'];

cd(dir_path)
addpath("Functions")

%%%%%%%%%%
%%%%%%%%%%

% Parameter structure
p = struct;

% Set parameters that are common uniform mobility (UM) & heterogenous
% mobility (HM) models
r = 0.04;                                  % discount rate 
p.ddelta = 0.0039;                         % exog. separations
p.bbeta = (1 - p.ddelta) * 1/(1+r)^(1/12); % discount factor       
p.eeta = 0.5;                              % elasticity of matching
p.z0 = 1;                                  % match starting productivity
p.F = 0;                                   % firing costs benchmark

%%%%

% Empirical targets
m = struct;
m.eu = 0.0173;               % EU rate
m.ue = 0.30;                 % UE ratee

%% 1. Uniform mobility (UM) model

% vectors of values for b
bb = linspace(0.55,0.95,50)';
L = length(bb);

% initialize vectors for outcomes
llambda = zeros(L,1);
eu = zeros(L,1);
ue = zeros(L,1);
urate = zeros(L,1);

% compute \lambda(b) and evaluate model over values for b
% using function "baseline_uniform"
for ii = 1:L

    p.b = bb(ii);
    [llambda(ii), el] =  model_UM(p,m);
    
    eu(ii) = el.eu;
    ue(ii) = el.ue;
    
end

% stack in a table
elasticities_UM = table(llambda,bb,eu,ue);

% display table
disp('Uniform-mobility (UM) elasticities:')
disp( elasticities_UM )

% save results
save( "Results.mat", "elasticities_UM" )


%% 2. Heterogenous mobility (HM) model

% initialize vectors for outcomes (including elasticity components)

% EU elasticity
eu_total = zeros(L,1);
eu_retention = zeros(L,1);
eu_composition = zeros(L,1);

% UE elasticity
ue_total = zeros(L,1);
ue_tightness = zeros(L,1);
ue_selection = zeros(L,1);

% Targets
eu = zeros(L,1);
eu_ratio = zeros(L,1);

% solve for values of lambda and sigma_x that min distance with targets

% initialize vectors
llambda = zeros(L,1);
sigma_x = zeros(L,1);

% starting point and bounds for algo
x0 = [0.10; 0.20];
lb = [0.05; 0.10];
ub = [0.20; 0.70];

% preallocate xstar
xstar = zeros(2,L);

% options for minimizer (fmincon, 'interior-point')
options = ...
optimoptions('fmincon', 'Display', 'final', 'Algorithm', 'interior-point', "ObjectiveLimit", 1.0000e-06, "MaxFunctionEvaluations", 500);

% functions to minimize
f_obj = cell(L,1);
p_ = cell(L,1);

for ii = 1:L
    
    p_{ii} = p;
    p_{ii}.b = bb(ii);
    f_obj{ii} =  @(x) model_HM(p_{ii}, false, x);

end

% loop over values for b
parfor ii = 1:L

% solve for parameters matching targets
xstar(:,ii) = fmincon(f_obj{ii}, x0, [], [], [], [], lb, ub, [], options);

% compute elasticities
[~, m, el] = model_HM(p_{ii}, true, xstar(:, ii));

eu_total(ii) = el.eu_total;
eu_retention(ii) = el.eu_retention;
eu_composition(ii) = el.eu_composition;

ue_total(ii) = el.ue_total;
ue_tightness(ii) = el.ue_tightness;
ue_selection(ii) = el.ue_selection

% moments
eu(ii) = m.eu;
eu_ratio(ii) = m.eu_ratio;

end

% gather parameters
for ii = 1:L
llambda(ii) = xstar(1,ii);
sigma_x(ii) = xstar(2,ii);
end

% stack in a table
elasticities_HM = table(bb, llambda, sigma_x, eu, eu_ratio, ...
      eu_retention, eu_composition, eu_total, ue_tightness, ue_selection, ue_total);

% display table
disp('Heterogeneous-mobility (HM) elasticities:')
disp( elasticities_HM )

% save results
save( "Results.mat", "elasticities_HM", '-append' )

%% End timer
total_time = toc;
disp('Baseline model analysis, total time:')
disp(toc)


