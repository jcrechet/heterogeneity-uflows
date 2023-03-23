function table_accounting = accounting_table(parameters, moments, data)


%% Initial structure
p0 = parameters;
m0 = moments;
clearvars parameters moments


%% Accounting for cc flows: contribution of policies and matching efficiency

% store some of the initial moments and paramaters
mean_W0 = m0.agg.mean_W;
mean_Y0 = m0.agg.mean_Y;

% initialize cells (6 entries: benchmark, policy, total, A, F, b + tax) 
J = 6;
p = cell(J,1);
m = cell(J,1);

% initialize counter
jj = 1;
p{jj} = p0;
m{jj} = m0;

% table with targets
targets = data.calibration_targets;

% targeted b_w ratio: equal to the ratio of b_w Eur vs US times the
% benchmark calibrated model's ratio
targets.b_w1 = 0.1 + m0.agg.b_w;

% targeted tax: pp difference Europe - US
targets.tax_w1 = 0.2;

% targeted firing costs
targets.f_w1 = targets.f_w_Eur;

% targeted ue and eu
targets.ue1 = ( targets.ue_Eur / targets.ue_US ) * m0.agg.ue;
targets.eu1 = ( targets.eu_Eur / targets.eu_US ) * m0.agg.eu;

% option for patternsearch
p0.pattern_options.MeshTolerance = 10^-3;

%% 1. Step 1: Calibrate firing costs, non-work utility, and tax to policy targets

step = 1;
disp('step 1: Calibrate firing costs, non-work utility, and tax to Europe policy targets')


%% step 1.1.: initial guess: estimates from partial equilibrium model

disp('step 1.1: assume partial equilibrium to search for initial guess for step 1.2 (calibration in gal eql)')
disp('')

% set options for calibration
p{jj}.equilibrium = false;
p{jj}.calibration = true;

% initial guess for the "initial-guess" search...
x0 = [ targets.f_w1 * mean_W0; ...
       targets.b_w1 * mean_W0; ...
       targets.tax_w1 * (mean_W0 / mean_Y0) ];

% bounds 
bd = [0.5*mean_W0, 1.5*mean_W0; 0, mean_W0; 0.10, 0.40];

% objective function
obj_fun = @(x) calibration_policies(x, p{jj}, targets, step);

% run minimization
x = patternsearch(obj_fun, x0, [], [], [], [], bd(:,1), bd(:,2), [], p0.pattern_options);

% update the initial guess
x0 = x;



%% step 1.2.: search in general equilibrium

disp(['step 1.2: general equilibrium search, using step 1.1. initial guess ' ...
    '(note: the search can hit some parameter points where no equilibrium is found!)'])
disp('')


% options
p{jj}.equilibrium = true;
p{jj}.display_equil_final = false;

% objective function
obj_fun = @(x) calibration_policies(x, p{jj}, targets, step);

% run minimization
x = patternsearch(obj_fun, x0, [], [], [], [], bd(:,1), bd(:,2), [], p0.pattern_options);

% update parameters
f1 = x(1);
b1 = x(2);
tax1 = x(3);

jj = jj + 1;
p{jj} = p{jj-1};
p{jj}.f = f1;
p{jj}.b = b1;
p{jj}.tax = tax1;
p{jj}.equilibrium = true;

% evaluate model at solution
disp('model outcomes after step 1:')
disp('')

[~, m{jj}] = equilibrium([], p{jj}, []);
disp(m{jj}.agg)


%% step 2. Determine matching efficiency that explain remaining US-Europe difference in the UE rate
step = 2;
disp('Step 2: Determine matching efficiency that explain remaining US-Europe difference in the UE rate')
disp('')

% objective function
obj_fun = @(x) calibration_policies(x, p{jj}, targets, step);

% step 2.1 rough grid search
L = 10;

% bounds for step 2.1 search
Alb = targets.ue1;
Aub = p{1}.A;
x = linspace(Alb, Aub, L)';

% search over grid
obj = zeros(L,1);
parfor ii = 1:L
    obj(ii) = obj_fun( x(ii) );
end

% bounds for next step
[~, i_min] = min( obj );
x0 = x(i_min);

% step 2.2: patternsearch
x = patternsearch(obj_fun, x0, [], [], [], [], Alb, Aub, [], p0.pattern_options);

% update parameters
A1 = x(1);
jj = jj+1;
p{jj} = p{jj-1};
p{jj}.A = A1;

% evaluate equilibrium
disp('final outcomes: ')
p{jj}.equilibrium = true;
[~, m{jj}] = equilibrium([], p{jj}, []);
disp(m{jj}.agg)
disp(' ')
disp('Next: decomposition to account for the contribution of different policies + matching efficiency')
disp('')

%% 3. Contribution of each factor, in isolation

% (i) firing costs
disp('Compute contribution of firing costs')
jj = jj + 1;
p{jj} = p{1}; p{jj}.f = f1;
p{jj}.equilibrium = true;
[~, m{jj}, ~] = equilibrium([], p{jj}, []);

% (ii) wedge from NWU and tax
disp('Compute contribution of NWU and taxes')
jj = jj + 1;
p{jj} = p{1}; p{jj}.b = b1; p{jj}.tax = tax1;
p{jj}.equilibrium = true;
[~, m{jj}, ~] = equilibrium([], p{jj}, []);

% (iv) matching efficiency
disp('Compute contribution of matching efficiency')
jj = jj + 1;
p{jj} = p{1}; p{jj}.A = A1;
p{jj}.equilibrium = true;
[~,  m{jj}, ~] = equilibrium([], p{jj}, []);


%% 4. Build table

disp('Build table with results')

table_accounting = table;

% labels
table_accounting.contribution = {'Benchmark'; 'Policy'; 'Total'; 'Firing costs'; 'NWU + Tax'; 'Matching efficiency'};


% initialize vectors

% flows
table_accounting.ue     =  zeros(J, 1);
table_accounting.eu     =  zeros(J, 1);

% parameter values
table_accounting.F      =  zeros(J, 1);
table_accounting.b      =  zeros(J, 1);
table_accounting.tax    =  zeros(J, 1);
table_accounting.A      =  zeros(J, 1);

% loop
for ii = 1:J
    
    % flows
    table_accounting.ue(ii)     = m{ii}.agg.ue;
    table_accounting.eu(ii)     = m{ii}.agg.eu;

    % parameters
    table_accounting.F(ii)      = p{ii}.f;
    table_accounting.b(ii)      = p{ii}.b;
    table_accounting.tax(ii)    = p{ii}.tax;
    table_accounting.A(ii)      = p{ii}.A;

end

end