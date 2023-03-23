%% Program for Europe-US differences
disp('Run programs for computing sources of secular unemployment Europe-US differences')

%% Load benchmark model parameters and statistics
load("Results\Benchmark.mat")
p = parameters; m = moments;
clearvars parameters moments

%% initialize seed state (for replication)
rng(1)

%% Calibration targets

% US u benefits replacement ratio
targets.b_w = m.agg.b_w;

% log-wage variance change
targets.var_w = m.agg.var_w * ( 1 /  ( 1 + data.calibration_targets.var_w_change ) );

% Europe policies
% Eur u benefits
targets.b_w1 = 0.1 + m.agg.b_w;

% Eur tax rate
targets.tax_w1 = 0.2;

% Eur firing costs
targets.f_w1 = 1;


%% Step 1: US labor market with low wage inequality/low complementarity

disp('Step 1: calibrate to US labor market, low wage inequality/lox complementarity (70s)')
disp(' ')

% 1.1. initial guess: estimates from partial equilibrium model
disp('Step 1.1.: gets an initial guess from partial equilibrium')
disp(' ')

p.equilibrium = false;

% initial points, bounds
x0 = [ p.b; p.chi ];
bd = [ 0, 2; 0.5, p.chi ];

% objective function
obj_fun = @(x) calibration_secular(x, p, targets);

% run minimization
p.pattern_options.MeshTolerance = 10^-3; 
x1 = patternsearch(obj_fun, x0, [], [], [], [], bd(:,1), bd(:,2), [], p.pattern_options);

% 1.2. step 2: general equilibrium
disp('Step 1.2.: general equilibrium')
disp(' ')

p.equilibrium = true;
p.display_equil_final = false;
obj_fun = @(x) calibration_secular(x, p, targets);
x2 = patternsearch(obj_fun, x1, [], [], [], [], bd(:,1), bd(:,2), [], p.pattern_options);

% update parameters using solution of step 1.2. 
p.b = x2(1);
p.chi = x2(2);

% compute new equilibrium
[~,m] = equilibrium([], p, []);

% save results
save('Results\alternative_models\secular.mat', "p", "m")


%% Step 2: European labor-market with low complementarity

disp('Step 2: calibrate to Europe labor market, low wage inequality/low complementarity (70s)')
disp(' ')

% load parameter structure from previous step
load('Results\alternative_models\secular.mat', 'p', 'm')

% 2.1. initial guess: estimates from partial equilibrium model
disp('Step 2.1.: partial equilibrium')
disp(' ')
p.equilibrium = false;
p.display_equil_final = false;

% initial guess and bounds
mean_W = m.agg.mean_W; mean_Y = m.agg.mean_Y; 
x0 = [ targets.f_w1 * mean_W; ...
       targets.b_w1 * mean_W; ...
       targets.tax_w1 * (mean_W / mean_Y) ];
bd = [0.5*mean_W, 3*mean_W; 0, mean_W; 0.10, 0.40];

% minimization
step = 1;
obj_fun = @(x) calibration_policies(x, p, targets, step);
x1 = patternsearch(obj_fun, x0, [], [], [], [], bd(:,1), bd(:,2), [], p.pattern_options);

% second step
disp('Step 2.2.: general equilibrium')
disp(' ')

p.equilibrium = true;
p.display_equil_final = false;
obj_fun = @(x) calibration_policies(x, p, targets, step);
x2 = patternsearch(obj_fun, x1, [], [], [], [], bd(:,1), bd(:,2), [], p.pattern_options);

% update parameters
p1 = p;
p1.f = x2(1);
p1.b = x2(2);
p1.tax = x2(3);

% compute new equilibrium
[~, m1] = equilibrium([], p1, []);

% save
save('Results\alternative_models\secular.mat', "p1", "m1", '-append')


%% Step 3: Impose European policies in the benchmark

disp('Step 3: Impose European policies in the benchmark US (90s-2000s) labor market')
disp(' ')

load('Results\Benchmark.mat')

p1 = parameters;
p1.f = table_accounting.F(2);
p1.b = table_accounting.b(2);
p1.tax = table_accounting.tax(2);

p1.equilibrium = true;
p1.display_equil_final = false;
[~, m1] = equilibrium([], p1, []);

% rename structures for convenience
p = parameters; m = moments; 

% save results
save('Results\Benchmark.mat', 'p', 'm', 'p1', "m1", '-append')


%% Table

table_secular = cell(2,1);
variable_names = {'US', 'Europe'};
row_names = {'\\chi', 'F', 'b', '\\phi', 'U rate', 'UE', 'EU', 'EE'};
table_size = [length(row_names) length(variable_names)];
vb_types = {'double', 'double'};

% benchmark
for jj = 1:2
    

    if jj == 1
        % benchmark
        load('Results\Benchmark.mat', 'p', 'm','p1', 'm1')

        % low complementarity environment
    else
        load('Results\alternative_models\secular.mat', 'p', 'm','p1', 'm1')
    end

table_secular{jj} = table( 'size', table_size, 'VariableTypes', vb_types, 'VariableNames', variable_names, 'RowNames', row_names);

table_secular{jj}{'\\chi', :}  = [ p.chi, p1.chi ]; 
table_secular{jj}{'F', :}      = [ p.f, p1.f ]; 
table_secular{jj}{'b', :}      = [ p.b, p1.b ]; 
table_secular{jj}{'\\phi', :}  = [ p.tax, p1.tax ]; 

table_secular{jj}{"U rate", :} = [m.agg.urate, m1.agg.urate]*100;
table_secular{jj}{"UE", :} = [m.agg.ue, m1.agg.ue]*100;
table_secular{jj}{"EU", :} = [m.agg.eu, m1.agg.eu]*100;
table_secular{jj}{"EE", :} = [m.agg.ee, m1.agg.ee]*100;

end


%% Latex table

t = table_secular;

fmt1 = '%.2g';
fmt2 = '%.2f';
line = @(ii, fmt)  [ ' &  &', num2str( t{1}{ii,1}, fmt ), ' & ',   num2str( t{1}{ii,2}, fmt ), ' &  &  &',   num2str( t{2}{ii,1}, fmt ), ' & ',  num2str( t{2}{ii,2}, fmt )  '\\\\' ];

Table_7 = [ ...
'\\begin{table}[!h]'        '\n' ...
'\\small'                   '\n' ...
'\\centering'               '\n' ...
'\\captionof{table}{Europe-U.S.\\ Secular unemployment differences}'                           '\n' ...
'\\label{tab:secular_urate}'                                                                   '\n' ...
'\\begin{tabular}{l ccc c ccc}'                                                                '\n' ...
'\\hline \\hline'                                                                              '\n' ...
'\\addlinespace'          '\n' ...
' & \\multicolumn{3}{c}{(i) \\textit{High complementarity} } &  & \\multicolumn{3}{c}{(ii) \\textit{Low complementarity} }   \\\\  ' '\n' ...
'\\addlinespace'          '\n' ...
'\\hspace{80pt}           & \\hspace{10pt}   & U.S. & Eur. &  & \\hspace{10pt}  & U.S. & Eur.             \\\\  '        '\n' ...
'\\addlinespace'          '\n' ...
'\\addlinespace'          '\n' ...
'\\textit{Parameter} &                        \\\\  '        '\n' ...
'\\addlinespace'          '\n' ...
'\\hspace{4pt}  $\\chi$ ',    line(1,fmt1), '\n' ...
'\\hspace{4pt}  $F$ ',        line(2,fmt1), '\n' ...
'\\hspace{4pt}  $b$ ',        line(3,fmt1), '\n' ...
'\\hspace{4pt}  $\\phi$ ',    line(4,fmt1), '\n' ...
'\\addlinespace'          '\n' ...
'\\addlinespace'          '\n' ...
'\\textit{Outcome} (\\%%)   &                 \\\\  '        '\n' ...
'\\addlinespace'          '\n' ...
'\\hspace{4pt}  U rate ',  line(5,fmt2), '\n' ...
'\\hspace{4pt}  UE ',      line(6,fmt2), '\n' ...
'\\hspace{4pt}  EU ',      line(7,fmt2), '\n' ...
'\\hspace{4pt}  EE ',      line(8,fmt2), '\n' ...
'\\addlinespace'          '\n' ...
'\\addlinespace'          '\n' ...
'\\hline \\hline'         '\n' ...
'\\end{tabular}'          '\n' ...
['\\caption*{\\footnotesize{\\textit{Note:} comparison of outcomes of the benchmark calibrated economy ' ...
'with high complementarity between skills and match quality as captured by $\\chi$ (i) and an economy with ' ...
'low complementarity (ii). Values of $\\chi$ are chosen to generate changes in log-wage variance consistent ' ...
'with \\cite{kambourov_manovskii_2009_Restud}. ' ...
'These economies are compared across different policy regimes: ' ...
'a ``U.S.'' regime (the benchmark calibrated policy parameters) and a ``Europe'' regime ' ...
'(see section 6.2. of the main text). The parameters are: $\\chi$: ' ...
'complementarity between skills and match quality; $F$: firing costs; $b$ non-work utility; $\\phi$: match output proportional tax. }}'] '\n' ...
'\\end{table}' ];

% Export table to latex files
tables_path = [results_path,'\tables\'];
fileID = fopen([tables_path, 'Table_7.tex'],'w');
fprintf(fileID, Table_7);
fclose(fileID);

%%

disp('Secular unemployment differences: done. Check latex files.')
