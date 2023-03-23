% Project: Heterogeneity in labor mobility and unemployment flows across countries
% Créchet, 2022
% Program name: "Alternative_models.m": 
% created November 2022; updated January 2023

% Compute and calibrate alternative models:
% model (i). no worker heterogeneity
% model (ii). no job heterogeneity
% model (iii). no worker/job heterogeneity
% model (i-b). no OJS
% model (iii-b). no worker/job het., no OJS

%% 1. Load benchmark model parameters and variables
load('Results/Benchmark.mat', "parameters", "moments")

tic

% benchmark paramaters
p0 = parameters;
clearvars parameters

% targets
targets = struct;
targets.ue = moments.agg.ue;
targets.eu = moments.agg.eu;
targets.ee = moments.agg.ee;

% targets = data.calibration_targets;

% benchmark replacement ratio
b_w0 = moments.agg.b_w;

% number of initial points for search
K = 20;

% bounds 
A_bd       = [0.2, 0.8]; 
b_bd       = [0.1, 3]; 
lambda_bd  = [0.01, 0.3];
s_bd       = [0.01, 1];

% matrix for bounds
bd = [A_bd; b_bd; lambda_bd; s_bd];

% initial points
tmp_rand = rand(length(bd), K);
x0 = tmp_rand.*bd(:,1) + (1-tmp_rand).*bd(:,2);

% OJS or not?
OJS = true;

% initialize seed state (for replication)
rng(1)

% option for patternsearch
p0.pattern_options.MeshTolerance = 10^-3;

%% Loop over alter models
disp('Compute and calibrate alternative models (i)-(iii)')

for jj = 1:3


    % Alternative model (i): no worker heterogeneity
    if jj == 1
        
        disp('Alternative model (i): no worker heterogneity')
        filename = 'model_1.mat';

        % initialize parameter structure
        p = p0;

        % set paramaters
        p.Ik = 1;
        p.Ix = 20;
        p.Iz = 10;


    % Alternative model (ii): no job heterogeneity
    elseif jj == 2

        disp('Alternative model (ii): no job heterogneity')
        filename = 'model_2.mat';

        p = p0;
        p.Ik = 10;
        p.Ix = 1;
        p.Iz = 20;


    % Alternative model (iii): no worker an no job heterogeneity
    elseif jj == 3

        disp('Alternative model (iii): no worker/no job heterogneity')
        filename = 'model_3.mat';

        p = p0;
        p.Ik = 1;
        p.Ix = 1;
        p.Iz = 20*10;

    end

% function to minimize
obj_fun = @(x) alt_model(x, p, targets, OJS);

% run minimization over grid of initial points
obj = zeros(K,1);
xx = 0*x0;

% pattern search over grid
for ii = 1:K
    disp('iteration:')
    disp(ii)
    [xx(:,ii), obj(ii)] = patternsearch(obj_fun, x0(:,ii), [], [], [], [], bd(:,1), bd(:,2), [], p0.pattern_options);
end

% select solutions
index_solutions = obj <= 1e-2;
solutions = xx(:,index_solutions);

% cell for evaluation in parallel
L = length(solutions(1,:));
p_tmp = cell(L,1);

for ii = 1:L
    p_tmp{ii} = p;
    p_tmp{ii}.equilibrium = false;
    p_tmp{ii}.calibration = true;
    p_tmp{ii}.A = solutions(1,ii);
    p_tmp{ii}.b = solutions(2,ii);
    p_tmp{ii}.lambda = solutions(3,ii);
    p_tmp{ii}.s = solutions(4,ii);
end

% vector to store b/Ew across solutions
b_w_tmp = zeros(1,L);

% evaluate solutions to compute b/Ew
parfor ii = 1:length( solutions(1,:) )

    [~, m] = equilibrium([], p_tmp{ii}, []);
    b_w_tmp(ii) = m.agg.b_w;

end

% identify solution with minimal distance with b/Ew in banchmark
[~, i_min] = min( abs( b_w_tmp - b_w0 ) );
xstar = solutions(:,i_min);

% update parameters
p.A = xstar(1);
p.b = xstar(2);
p.lambda = xstar(3);
p.s = xstar(4);

% evaluate model
p.equilibrium = false;
[~, m, v, ~] = equilibrium([], p, []);
disp(m.agg)

% check equilibrium
p.c_v = v.c_v;
p.equilibrium = true;
equilibrium([], p, []);

% save results
parameters = p;
moments = m;
variables = v;

save(['Results\alternative_models\', filename], "parameters", "variables", "moments" )

clearvars p m v i_min b_w_tmp p_tmp xx obj

end


%% alternative models: no OJS 
disp('Compute and calibrate alternative models (i-b) and (iii-b)')

OJS = false;
p0.s = 0;
   
% matrix for bounds
A_bd = [0.3, 0.8];
b_bd = [0, 1];
lambda_bd = [0.01, 0.25];
bd = [A_bd; b_bd; lambda_bd];

% initial points
tmp_rand = rand(length(bd), K);
x0 = tmp_rand.*bd(:,1) + (1-tmp_rand).*bd(:,2);

for jj = 1:2


    % Alternative model (i-b): no OJS
    if jj == 1

        disp('Alternative model (i-b): no OJS')
        filename = 'model_1b.mat';

        % initialize parameter structure
        p = p0;

    % Alternative model (iii-b): no OJS and no heterogeneity
    elseif jj == 2

        disp('Alternative model (iii-b): no heterogeneity, no OJS')
        filename = 'model_3b.mat';

        p = p0;
        p.Ik = 1;
        p.Ix = 1;
        p.Iz = 20*10;

    end

    % function to minimize
    obj_fun = @(x) alt_model(x, p, targets, OJS);

    % run minimization over grid of initial points
    obj = zeros(K,1);
    xx = 0*x0;

    % pattern search over grid
    for ii = 1:K
        [xx(:,ii), obj(ii)] = patternsearch(obj_fun, x0(:,ii), [], [], [], [], bd(:,1), bd(:,2), [], p0.pattern_options);
    end

    % select candidate solutions
    index_solutions = obj < 1e-3;
    solutions = xx(:, index_solutions);

    % cell for evaluation in parallel
    L = length(solutions(1,:));
    p_tmp = cell(L,1);

    for ii = 1:L
        p_tmp{ii} = p;
        p_tmp{ii}.equilibrium = false;
        p_tmp{ii}.calibration = true;
        p_tmp{ii}.A = solutions(1,ii);
        p_tmp{ii}.b = solutions(2,ii);
        p_tmp{ii}.lambda = solutions(3,ii);
    end


    % vector to store b/Ew across solutions
    b_w_tmp = zeros(1,L);

    % evaluate solutions to compute b/Ew
    parfor ii = 1:length( solutions(1,:) )

        [~, m] = equilibrium([], p_tmp{ii}, []);
        b_w_tmp(ii) = m.agg.b_w;

    end

    % identify solution with minimal distance with b/Ew in banchmark
    [~, i_min] = min( abs( b_w_tmp - b_w0 ) );
    xstar = solutions(:,i_min);

    % update parameters
    p.A = xstar(1);
    p.b = xstar(2);
    p.lambda = xstar(3);

    % evaluate model
    p.equilibrium = false;
    [~, m, v, ~] = equilibrium([], p, []);
    disp(m.agg)

    % check equilibrium
    p.c_v = v.c_v;
    p.equilibrium = true;
    equilibrium([], p, []);

    % save results
    parameters = p;
    moments = m;
    variables = v;

    save(['Results\alternative_models\', filename], "parameters", "variables", "moments" )

    clearvars p m v i_min b_w_tmp p_tmp xx obj

end

%%

disp('Alternative_model, time:')
disp(toc)

