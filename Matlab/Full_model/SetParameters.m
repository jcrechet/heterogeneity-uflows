% This script creates a structure containing parameters that are passed as an argument to
% functions solving numerically the model's equilibrium and conducting the
% numerical minimization in the calibration step 

%% 

% initialize structure called "p" 
p = struct;

% structure for timer
p.timer = [];


%% 1. Settings for numerical solution, simulations, and quantitative experiments

% logical variables: solve for partial (for calibration/validation) or general
% equilibrium (for experiments)
p.equilibrium = false;

% set true to display equilibrium computation progress in all iterations
p.display_equil_iterations = false;
p.display_equil_final = true;

% set true to compute only moments used in calibration
p.calibration = true; 

% set true to compute steady-state equilibrium distribution of state
% variables (for computing equilibrium or the vacancy posting cost consistent with theta = 1 and calibrated parameters)
p.sstate = false;

% max nb of iterations on equil tightness in fixed-point search
p.max_iter_eq = 100;

% criterion (maximum distance) for stopping search of a fixed-point for equilibrium tightness
p.crit_eq = 1.0000e-05;

% logical variables for alternative models (if "false": shut down the
% associated channel)
p.skills = true;               % skill heterogeneity
p.matchquality = true;         % match quality
p.ojs = true;                  % OJS
p.searchInt = true;            % variable search intensity 

% time unit
p.periodlen = 1/12; % 1/12: monthly; 1/26: bi-weekly, etc.

% entry age in years
p.age_entry = 19;

% duration of working life in years
p.working_life = 43;

% retirement age in model time units
p.T = p.working_life*(1/p.periodlen);

% age range for computing aggregate stats
p.age_range = [20, 60];

% grid sizes
p.Iz = 5;  % stochastic match productivity
p.Ix = 11; % permanent match quality (pick odd number)
p.Ik = 7;  % ability and human capital grid


%% 2. Model preset parameters

% 2.1. Preset parameters

p.Beta = (1.04)^(-p.periodlen); % Discount factor

p.Gamma = 0.3;                  % Bargaining power

p.Delta = 0.0039;               % Exog. separation rates

p.Eta = 0.5;                    % Elasticity of meeting

p.chi = 1;                      % Complementarity bw worker skills and match quality

p.c1 = 1;                       % Curvature of search function  

p.f = 0;                        % Firing costs

p.tax = 0;                      % Proportional tax on output

p.c_v = 1.683179402135146;      % Vacancy posting cost
  
p.mu0 = 0;                      % Expected skill growth rate in unemployment                      


%% 3. Model internally calibrated parameters (by SMM; these are the benchmark/initial values)

% This part specifies initial values for the parameters that are internally
% calibrated, and stack them in a vector that is passed to the numerical
% model and to the objective function upon which the calibration
% minimization algorithm relies. It also specifies a 2-col matrix with bounds associated with each parameter. 
% Note that each parameter initial value and bounds can be recovered using its
% associated index value. This indexing scheme is also used in the function
% "numerical model" to associate each model parameter to its appropriate value in the parameter vector. 

% initialize index referring to the vector components
j=1;

% non-work income
p.b =  1.371435356752773;          
p.x0(j) = p.b;      p.ib=j;                 p.bd(j,:) = [0.1,2];        p.label{j} = 'b';        j=j+1;

% matching efficiency
p.A = 0.586807770793985;        
p.x0(j) = p.A;     p.iA=j;                  p.bd(j,:) = [0.4,0.8];      p.label{j} = 'A';        j=j+1;

% linear component of search cost function
p.c0 = 1.717513273754121;
p.x0(j) = p.c0;   p.ic0=j;                  p.bd(j,:) = [0.5,10];        p.label{j} = '\\chi_e';   j=j+1;

% on-the-job search efficiency
p.s = 0.564024813298746;  
p.x0(j) = p.s;  p.is=j;                     p.bd(j,:) = [0.1,1];      p.label{j} = 's';           j=j+1;

% probability of stochastic shock
p.lambda = 0.193243493568808;
p.x0(j) = p.lambda; p.ilambda=j;            p.bd(j,:) = [0.05,0.30];    p.label{j} = '\\lambda';    j=j+1;

% std of ability disturbance
p.sigma_a = 0.089180017324280;
p.x0(j) = p.sigma_a; p.isigma_a=j;          p.bd(j,:) = [0.01,0.30];    p.label{j} = '\\lambda';    j=j+1;

% persistence of ability 
p.rho_a = 0.966718153711846;
p.x0(j) = p.rho_a; p.irho_a=j;              p.bd(j,:) = [0.95,1-10^(-4)];  p.label{j} = '\\lambda';    j=j+1;

% std of log permanent match quality
p.sigma_x = 0.446392846966210;     
p.x0(j) = p.sigma_x; p.isigma_x=j;          p.bd(j,:) = [0.20,0.60];       p.label{j} = '\\sigma_y';   j=j+1;

% skill growth rate, employment
p.mu1 = 0.064027503384152;
p.x0(j) = p.mu1; p.imu1=j;                  p.bd(j,:) = [0.01,0.1];        p.label{j} = '\\mu_1';       j=j+1;  

% free search effort
p.slb = 0.540427803807512;
p.x0(j) = p.slb;    p.islb=j;               p.bd(j,:) = [0.0,0.6];         p.label{j} = '\\underline{s}';  %j=j+1;

% Total number of parameters
p.K = j;

p.x0 = p.x0';
p.lb = p.bd(:,1); p.ub = p.bd(:,2);
p.label = p.label';

%% Pattern search options for calibration
timelimit = 2*60*60; 
p.pattern_options = optimoptions('patternsearch', 'PlotInterval', 10, 'PlotFcns', {@psplotbestx, @psplotfuncount, @psplotmeshsize, @psplotmeshsize}, ...
    'Display', 'iter', 'UseParallel', 'always', 'TimeLimit', timelimit, 'MeshTolerance', 10^-3 );

%% rename structure
parameters = p;
clearvars p


%% cleanup
clearvars j timelimit
