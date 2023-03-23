% Project: Heterogeneity in labor mobility and unemployment flows across countries
% Créchet, 2022 (revised 2023)
% Program name: "Elasticities.m"
% created January 2023; last updated March 2023

% Compute decomposition of elasticities and variation of outcomes (urate
% and aggregat flows) with respect to policies and matching efficiency


%% 0. load benchmark model parameters and varibles
load('Results/Benchmark.mat')

% rename benchmark object
p0 = parameters;
m0 = moments;
v0 = variables;

clearvars parameters moments variables

% mean salary equilibrium benchmark
mean_W0 = m0.agg.mean_W;

% benchmark parameters
b0 = p0.b;
A0 = p0.A;

% select computation options
p0.equilibrium = true;
p0.calibration = false;

% preset parameters 
L = 20; % lenght of vector for counterfactual param values
K = p0.K; % lenght of vector for model calibrated parameters


%% 1. Compute counterfactual model outcomes for different policy values (in isolation)

% initialize vectors for counterfactual policies 

% 1.1. Firing costs
f = linspace(0, 3, L)';     % vector for firing costs in terms of benchmar salary

% 1.2. Non-work utility
b = linspace(1, 1.15, L)';  % in terms of benchmark param value

% 1.3. output tax
tax_y = linspace(0, 0.15, L)'; % percentage of output

% 1.4. Efficiency of matching
A = linspace(1, 0.7, L)';   % expressed in terms of benchmark param value

% Loop over types of policy
for jj = 1:4


    % initialize structures
    % parameters
    p = cell(L,1);
    for ii = 1:L
        p{ii} = p0;
    end

    % equilibrium variables
    v = cell(L,1);
    v{1} = v0;

    % model moments
    m = cell(L,1);
    m{1} = m0;

    % vector of mean wage
    mean_W = zeros(L,1);
    mean_W(1) = mean_W0;


    
    % Set policy-specific variables
    
    % 1. Firing costs
    if jj == 1
        
        % update param structure
        for ii = 2:L
            p{ii}.f = f(ii)*mean_W0;
        end
        
    % 2. non-work utility
    elseif jj == 2

        for ii = 2:L
            p{ii}.b = b(ii)*b0;
        end

    % 3. non-work income
    elseif jj == 3
        
        for ii = 2:L
            p{ii}.tax = tax_y(ii);
        end

    % 4. Efficiency of matching
    else

        for ii = 2:L
            p{ii}.A = A(ii)*A0;
        end

    end


    % Compute counterfactual outcomes
    parfor ii = 2:L
        [~, m{ii}, v{ii} ] = equilibrium([], p{ii}, []);
        mean_W(ii) = m{ii}.agg.mean_W;
    end



    % compute implied policy vector in terms of mean wage (if applicable)
    % and save outcomes

    % Firing costs
    if jj == 1
        
        f_w = (f*mean_W0) ./ mean_W; 
        save('Results\decompositions\firing_costs.mat', "v", "m", "p", "f", "f_w")

    % Non-work income
    elseif jj == 2
        
        b_w = (b*b0) ./ mean_W;
        save('Results\decompositions\non_work_utility.mat', "v", "m", "p", "b", "b_w")
    
    % tax
    elseif jj == 3
        
        mean_Y = zeros(L,1);
        for ii = 1:L
            mean_Y(ii) = m{ii}.agg.mean_Y;
        end
        tax_w = ( tax_y .* mean_Y )  ./ mean_W;

        save('Results\decompositions\output_tax.mat', "v", "m", "p", "tax_y", "tax_w")

    % matching efficiency
    else

        save('Results\decompositions\matching_efficiency.mat', "v", "m", "p", "A")

    end

end

% cleanup
clearvars p0 m0 v0 L  K f b tax_y A mean_W mean_y f_w b_w b0 A0


%% 2. Compute variations in aggregate rates and their components

% loop over policies
for jj = 1:4

    % 2.1. firing costs
    if jj == 1
    load('Results\decompositions\firing_costs.mat', "v", "m", "f", "f_w", "p")

    % 2.2. non-work utility
    elseif jj == 2
    load('Results\decompositions\non_work_utility.mat', "v", "m", "b", "b_w", "p")

    % 2.3. tax output
    elseif jj == 3
    load('Results\decompositions\output_tax.mat', "v", "m", "tax_y", "tax_w", "p")

    % 2.4. matching efficiency
    elseif jj == 4
    load('Results\decompositions\matching_efficiency.mat', "v", "m", "A", "p")
    end


% preallocate tables
L = length(m);

% ue elasticity and components
ue = array2table( zeros(L,7), "VariableNames", ...
    ["actual", "total", "residual", "skill", "tightness", "search", "selection"]);

% eu 
eu = array2table( zeros(L,7), "VariableNames", ...
    ["actual", "total", "residual", "skill", "job", "retention", "reallocation"]);

% ee
ee = array2table( zeros(L,7), "VariableNames", ...
    ["actual", "total", "residual", "skill", "job", "contact", "selection"]);


% a) UE rate decomposition
for ii = 2:L
ue(ii,:) = ue_variation(m{1}, m{ii}, v{1}, v{ii}, p{1});
end

% b) EU rate decomposition
for ii = 2:L
eu(ii,:) = eu_variation(m{1}, m{ii}, v{1}, v{ii}, p{1});
end

% c) EE rate
for ii = 2:L
ee(ii,:) = ee_variation(m{1}, m{ii}, v{1}, v{ii}, p{1});
end

% d) table with aggregate moments
agg = array2table( zeros(L,4), "VariableNames", ...
    ["u", "ue", "eu", "ee"]);
for ii = 1:L
    agg.u(ii) = m{ii}.agg.urate;
    agg.ue(ii) = m{ii}.agg.ue;
    agg.eu(ii) = m{ii}.agg.eu;
    agg.ee(ii) = m{ii}.agg.ee;
end


% add variables for policy values and save

% firing costs
if jj == 1
    
    ue = addvars(ue, f, f_w, 'Before', 'actual');
    eu = addvars(eu, f, f_w, 'Before', 'actual');
    ee = addvars(ee, f, f_w, 'Before', 'actual');

    save('Results\decompositions\firing_costs.mat', "ue", "eu", "ee", "agg") 

% non-work utility
elseif jj == 2

    ue = addvars(ue, b, b_w, 'Before', 'actual'); 
    eu = addvars(eu, b, b_w, 'Before', 'actual');
    ee = addvars(ee, b, b_w, 'Before', 'actual');

    save('Results\decompositions\non_work_utility.mat', "ue", "eu", "ee", "agg") 


% output tax
elseif jj == 3

    ue = addvars(ue, tax_y, tax_w, 'Before', 'actual');
    eu = addvars(eu, tax_y, tax_w, 'Before', 'actual');
    ee = addvars(ee, tax_y, tax_w, 'Before', 'actual');

    save('Results\decompositions\output_tax.mat', "ue", "eu", "ee", "agg") 


% matching efficiency
elseif jj == 4

    ue = addvars(ue, A, 'Before', 'actual');
    eu = addvars(eu, A, 'Before', 'actual');
    ee = addvars(ee, A, 'Before', 'actual');

    save('Results\decompositions\matching_efficiency.mat', "ue", "eu", "ee", "agg") 

end

end


%% 4. Unemployment change decomposition

% load results
for jj = 1:4

    if jj == 1 
    load('Results\decompositions\firing_costs.mat', "ue", "eu", "agg")
    elseif jj == 2
    load('Results\decompositions\non_work_utility.mat', "ue", "eu", "agg") 
    elseif jj == 3
    load('Results\decompositions\output_tax.mat', "ue", "eu", "agg") 
    elseif jj == 4
    load('Results\decompositions\matching_efficiency.mat', "ue", "eu", "agg") 
    end


% initialize table
u = table;

% unemployment change
u.actual = ( agg.u - agg.u(1) ) ./ agg.u(1);

% contribution of UE rate change
u.ue = (1/2)*( agg.eu ./ (agg.eu + agg.ue) -  agg.eu ./ (agg.eu + agg.ue(1)) ...
               + agg.eu(1) ./ (agg.eu(1) + agg.ue) - agg.eu(1) ./ (agg.eu(1) + agg.ue(1)) ) ./ agg.u(1);

% contribution of EU rate change
u.eu = (1/2)*( agg.eu ./ (agg.eu + agg.ue) -  agg.eu(1) ./ (agg.eu(1) + agg.ue) ...
               + agg.eu ./ (agg.eu + agg.ue(1)) - agg.eu(1) ./ (agg.eu(1) + agg.ue(1)) ) ./ agg.u(1) ;


% contribution of distribution
u.distr_ue = ( ue.skill ./ ue.total ) .* u.ue;              % throught UE
u.distr_ue(1) = 0;

u.distr_eu = ( (eu.skill + eu.job) ./ eu.total ) .* u.eu;   % throught EU 
u.distr_eu(1) = 0;

u.distr = u.distr_ue + u.distr_eu;                          % total 


% contribution of conditional transition

% throught UE
tmp = ue.tightness + ue.search + ue.selection;
u.trans_ue = ( tmp ./ ue.total ) .* u.ue;
u.trans_ue(1) = 0;

% throught EU
tmp = eu.retention + eu.reallocation;
u.trans_eu = ( tmp ./ eu.total ) .* u.eu;
u.trans_eu(1) = 0;

% total 
u.trans = u.trans_ue + u.trans_eu;


% total computed
u.total = u.distr + u.trans;
u.residual = u.actual - u.total; 


% append table, save
if jj == 1

    f = ue.f; f_w = ue.f_w;
    u = addvars( u, f, f_w, 'Before', 'actual' );
    save('Results\decompositions\firing_costs.mat', "u", '-append')

elseif jj == 2

    b = ue.b; b_w = ue.b_w;
    u = addvars( u, b, b_w, 'Before', 'actual' );
    save('Results\decompositions\non_work_utility.mat', "u", '-append')

elseif jj == 3

    tax_y = ue.tax_y; tax_w = ue.tax_w;
    u = addvars( u, tax_y, tax_w, 'Before', 'actual' );
    save('Results\decompositions\output_tax.mat', "u", '-append')

elseif jj == 4

    A = ue.A;
    u = addvars( u, A, 'Before', 'actual' );
    save('Results\decompositions\matching_efficiency.mat', "u", '-append')

end

end


%% 4. Elasticity tables


% initialize matrix for values to fill tables
pp = zeros(6*3+9,4);

for jj = 1:4

    % firing costs
    if jj == 1

        % load results
        load('Results\decompositions\firing_costs.mat', "ue", "eu", "ee", "u")

        % variation in firing costs: one unit of the numéraire (same as
        % baseline model)
        tmp = ue.f;
        [~, kk] = min( abs( tmp - 1 ) );
        var_policy = tmp(kk) / 100;
        clearvars tmp

    % non-work utility
    elseif jj == 2
        
        % results
        load('Results\decompositions\non_work_utility.mat', "ue", "eu", "ee", "u")
        
        % variation: one p.p. increase in replacement ratio
        tmp = ue.b_w - ue.b_w(1);
        [~, kk] = min( abs( tmp - 0.01 ) ); kk = max( kk, 2 );
        var_policy = tmp(kk);
        clearvars tmp

    % proportional tax
    elseif jj == 3

        % results
        load('Results\decompositions\output_tax.mat', "ue", "eu", "ee", "u")

        % variation: one p.p. increase in proportional wage tax
        tmp = ue.tax_w;
        [~, kk] = min( abs( tmp - 0.01 ) ); kk = max( kk, 2 );
        var_policy = tmp(kk);
        clearvars tmp

    % efficiency of matching
    elseif jj == 4

        % results
        load('Results\decompositions\matching_efficiency.mat', "ue", "eu", "ee", "u")

        % variation: one percentage point reduction in matching efficiency
        tmp = ( ue.A - ue.A(1) ) / ue.A(1);
        [~, kk] = min( abs( tmp - 0.01 ) ); kk = max( kk, 2 );
        var_policy = - tmp(kk);
        clearvars tmp

    end

% ue rate elasticity
pp(1:6,jj) = [ ue.skill(kk); ...
               ue.tightness(kk) + ue.search(kk) + ue.selection(kk); ...
               ue.tightness(kk); ue.search(kk); ue.selection(kk); ...
               ue.total(kk) ] / var_policy;

% eu rate
pp(7:12,jj) = [eu.skill(kk); eu.job(kk); ...
               eu.retention(kk) + eu.reallocation(kk); ...
               eu.retention(kk); eu.reallocation(kk); ...
               eu.total(kk) ] / var_policy;

% ee rate
pp(13:18,jj) = [ee.skill(kk); ee.job(kk); ...
                ee.contact(kk) + ee.selection(kk); ...
                ee.contact(kk); ee.selection(kk); ...
                ee.total(kk) ] / var_policy;

% u rate
pp(19:end,jj) = [u.ue(kk); u.distr_ue(kk); u.trans_ue(kk); ...
                 u.eu(kk); u.distr_eu(kk); u.trans_eu(kk); ...
                 u.distr(kk); u.trans(kk); u.total(kk) ] / var_policy;

end

% set formats for entries in table
format1 = '%.1f';
format2 = '%.1f';

% function to fill in table lines with data from matrix with results
% plain
line = @(ll) [ ' ' num2str( pp(ll,1), format1), ' & ', num2str( pp(ll,2), format1), ' & ', num2str( pp(ll,3), format1), ' & ', num2str( pp(ll,4), format1), ' \\\\'];

% italic
line_it = @(ll) [ ' \\textit{', num2str( pp(ll,1), format1), '} & ', ' \\textit{',num2str( pp(ll,2), format1), '} & ', ' \\textit{',num2str( pp(ll,3), format1), '} & ', ' \\textit{',num2str( pp(ll,4), format1),'} ', '\\\\'];

% 4.1. Aggregate flows

Table_3 =  ...
    [ '\\begin{table}[!h]' '\n' ...
    '\\centering' '\n' ...
    '\\caption{Decomposition of quantitative model steady-state flow semi-elasticities}' '\n' ...
    '\\label{tab:table_3}'                                         '\n' ...
    '\\small'                                                                       '\n' ...
    '\\begin{tabular}{l l c c c c}'                                                 '\n' ...
    '\\hline \\hline'                                                               '\n' ...
    '\\addlinespace'                                                                '\n' ...
    '  & \\hspace{80pt}  & \\hspace{12.5pt} $F$ \\hspace{12.5pt} & \\hspace{12.5pt} $b$ \\hspace{12.5pt} & \\hspace{12.5pt} $\\phi$ \\hspace{12.5pt} & \\hspace{12.5pt} $A$ \\hspace{12.5pt} \\\\'  '\n' ...
    '\\addlinespace'  '\n' ...
    '\\addlinespace'  '\n' ...
    '\\text{UE elasticity contribution (\\%%)}' '\\\\'                                                      '\n' ...
    '\\addlinespace'                                                              '\n' ...
    '\\hspace{2pt}     skill distribution           &   &',        line(1),       '\n' ...
    '\\hspace{2pt}     conditional transition       &   &',        line(2),       '\n' ...
    '\\hspace{6pt}     \\textit{tightness}          &   &',        line_it(3),    '\n' ...
    '\\hspace{6pt}     \\textit{search}             &   &',        line_it(4),    '\n' ...
    '\\hspace{6pt}     \\textit{selection}          &   &',        line_it(5),    '\n' ...
    '\\addlinespace'                                                              '\n' ...
    'Total                                          &   &',        line(6)        '\n' ...
    '\\addlinespace'                                                              '\n' ...
    '\\addlinespace'                                                              '\n' ...
    '\\addlinespace'                                                              '\n' ...
    '\\text{EU elasticity contribution (\\%%)}' '\\\\'                                                      '\n' ...
    '\\addlinespace'                                                              '\n' ...
    '\\hspace{2pt}     skill distribution           &   &'         line(7)        '\n' ...
    '\\hspace{2pt}     job distribution             &   &'         line(8)        '\n' ...
    '\\hspace{2pt}     conditional transition       &   &'         line(9)        '\n' ...
    '\\hspace{6pt}     \\textit{retention}          &   &'         line_it(10)    '\n' ...
    '\\hspace{6pt}     \\textit{(no) reallocation}  &   &'         line_it(11)    '\n' ...
    '\\addlinespace'                                                              '\n' ...
    'Total                                          &   &'         line(12)       '\n' ...
    '\\addlinespace'                                                              '\n' ...
    '\\addlinespace'                                                              '\n' ...
    '\\addlinespace'                                                              '\n' ...
    '\\text{EE elasticity contribution (\\%%)}' '\\\\'                            '\n' ...
    '\\addlinespace'                                                              '\n' ...
    '\\hspace{2pt}     skill distribution           &   &'        line(13)        '\n' ...
    '\\hspace{2pt}     job distribution             &   &'        line(14)        '\n' ...
    '\\hspace{2pt}     conditional transition       &   &'        line(15)        '\n' ...
    '\\hspace{6pt}     \\textit{contact}            &   &'        line_it(16)     '\n' ...
    '\\hspace{6pt}     \\textit{selection}          &   &'        line_it(17)     '\n' ...
    '\\addlinespace'                                                              '\n' ...
    'Total                                          &   &'        line(18)        '\n' ...
    '\\addlinespace'                                                              '\n' ...
    '\\addlinespace'                                                              '\n' ...
    '\\addlinespace'                                                              '\n' ...
    '\\hline \\hline'                                                             '\n' ...
    '\\end{tabular}'         '\n' ...
    ['\\caption*{\\footnotesize{\\textit{Notes:} decomposition of quantitative model aggregate steady-state UE, EU, and EE ' ...
    'semi-elasticities with respect to the following parameters: $F$: firing costs; $b$ non-work utility; $\\phi$: match-output proportional tax; ' ...
    '$A$: matching efficiency. The decomposition is based on \\eqref{model:UE_elasticity}, ' ...
    '\\eqref{model:EU_elasticity} and \\eqref{appendix:EE_elasticity}. Each entry represents the percentage-point value of a component in the total elasticity. ' ...
    'The semi-elasticities take as reference: a unit numéraire increase in $F$; a one-percentage point increase in the replacement ratio $\\zeta = b / E(w)$ and the effective tax rate $\\phi_w = \\phi E(f(x,y,z)) / E(w)$; a one percent decrease in matching efficiency.' ...
    '}}'],  '\n' ...
    '\\end{table}' ];


% Export table
tables_path = [results_path,'\Tables\'];
fileID = fopen([tables_path, 'Table_3.tex'],'w');
fprintf(fileID, Table_3);
fclose(fileID);


% 4.2. Unemployment rate
Table_4 =  ...
    [ '\\begin{table}[!h]' '\n' ...
    '\\centering' '\n' ...
    '\\caption{Decomposition of quantitative model steady-state unemployment semi-elasticities}' '\n' ...
    '\\label{tab:table_4}'                       '\n' ...
    '\\small'                                                                               '\n' ...
    '\\begin{tabular}{l l c c c c}'                                                         '\n' ...
    '\\hline \\hline'                                                                        '\n' ...
    '\\addlinespace'                                                                         '\n' ...
    ' & \\hspace{80pt}  & \\hspace{12.5pt} $F$ \\hspace{12.5pt}  & \\hspace{12.5pt} $b$  \\hspace{12.5pt} & \\hspace{12.5pt} $\\phi$ \\hspace{12.5pt} & \\hspace{12.5pt} $A$ \\hspace{12.5pt} \\\\'  '\n' ...
    '\\addlinespace'                                                                         '\n' ...
    '\\addlinespace'                                                                         '\n' ...
    'UE contribution (\\%%)                               &   &',       line(19),       '\n' ...
    '\\hspace{4pt}     \\textit{distribution}             &   &',       line_it(20),    '\n' ...
    '\\hspace{4pt}     \\textit{conditional transition}   &   &',       line_it(21),    '\n' ...
    '\\addlinespace'  '\n' ...
    'EU contribution (\\%%)                               &   &',       line(22),       '\n' ...
    '\\hspace{4pt}     \\textit{distribution}             &   &',       line_it(23),    '\n' ...
    '\\hspace{4pt}     \\textit{conditional transition}   &   &',       line_it(24),    '\n' ...
    '\\addlinespace'  '\n' ...
    'Distribution, total                                  &   &',       line(25),       '\n' ...
    'Transition, total                                    &   &',       line(26),       '\n' ...
    '\\addlinespace'  '\n' ...
    'Total unemployment elasticity                        &   &',       line(27),       '\n' ...
    '\\addlinespace' '\n' ...
    '\\hline \\hline' '\n' ...
    '\\end{tabular}'         '\n' ...
    ['\\caption*{\\footnotesize{\\textit{Notes:} decomposition of quantitative model aggregate steady-state unemployment ' ...
    'semi-elasticities with respect to the following parameters: $F$: firing costs; $b$ non-work utility; $\\phi$: match-output proportional tax; ' ...
    '$A$: matching efficiency. The decomposition is based on \\eqref{model:UE_elasticity}, ' ...
    '\\eqref{model:EU_elasticity} and \\eqref{appendix:urate_elasticity}. Each entry represents the percentage-point value of a component in the total elasticity. ' ...
    '``Distribution, total'''', and ``transition, total'''' sum up the components associated with the UE and EU elasticities. ' ...
    'The semi-elasticities take as reference: a unit numéraire increase in $F$; a one-percentage point increase in the replacement ratio $\\zeta = b / E(w)$ and the effective tax rate $\\phi_w = \\phi E( f(x,y,z) / E(w)$; a one percent decrease in matching efficiency.' ...    
    '}}' ...
    ],  '\n' ...
    '\\end{table}' ];

% Export table
fileID = fopen([tables_path, 'Table_4.tex'],'w');
fprintf(fileID, Table_4);
fclose(fileID);


%% 5. Figures

for jj = 1:4


    % firing costs
    if jj == 1

        load('Results\decompositions\firing_costs.mat', "ue", "eu", "ee", "u")
        indepvar = ue.f_w;
        filename = 'Figure_5';
        x_label = 'Firing costs ($F$) / average equil. wage (E(w))';

        % non-work utility
    elseif jj == 2

        load('Results\decompositions\non_work_utility.mat', "ue", "eu", "ee", "u")
        indepvar = ue.b;
        filename = 'Figure_6';
        x_label = 'Non-work utility rel.\ to benchmark';

        % tax
    elseif jj == 3

        load('Results\decompositions\output_tax.mat', "ue", "eu", "ee", "u")
        indepvar = ue.tax_w;
        filename = 'Figure_7';
        x_label = 'Equilibrium average wage tax rate $\phi_w$';

        % matching efficiency
    elseif jj == 4

        load('Results\decompositions\matching_efficiency.mat', "ue", "eu", "ee", "u")
        indepvar = ue.A;
        filename = 'Figure_8';
        x_label = 'Matching efficiency ($A$) rel.\ to benchmark';

    end

% UE rate
fig = figure;

fig.Position = [1000 250 1000 600]; 

subplot(2,2,1)

hold on

% total
plot(indepvar, ue.total, 'Color', '#0D16A0', 'LineWidth', 2.25);
% skill distribution
plot(indepvar, ue.skill, 'Color', '#5AA2E3', 'LineStyle', '-.', 'LineWidth', 1.5);
% contact
plot(indepvar, ue.tightness + ue.search, 'Color', '#5AA2E3', 'LineStyle', ':', 'LineWidth', 1.5);
% selection
plot(indepvar, ue.selection, 'Color', '#5AA2E3', 'LineStyle', '--', 'LineWidth', 1.5);
% zero line
yline(0, ':k', 'LineWidth', 1.00)

xlim( [min(indepvar) max(indepvar)] )
grid("on") 
box('on')

title('(a) UE rate', 'Interpreter', 'latex', 'Interpreter', 'latex')
legend('total', 'skill distribution', 'contact', 'selection', '','Interpreter','latex');
legend('boxoff')

if jj == 1
    ylim([ -0.8 0.6 ])
    yticks([-0.75 -0.50 -0.25 0 0.25 0.5])
    legend('location', 'northwest')
else
    legend('location', 'best')
end
hold off

% EU rate
subplot(2,2,2)
hold on

% total
plot(indepvar, eu.total, 'Color', '#0D16A0', 'LineWidth', 2.25);
% distribution
plot(indepvar, eu.skill + eu.job, 'Color', '#5AA2E3', 'LineStyle', '-.', 'LineWidth', 1.5);
% retention
plot(indepvar, eu.retention, 'Color', '#5AA2E3', 'LineStyle', ':', 'LineWidth', 1.5);
% realocation
plot(indepvar, eu.reallocation, 'Color', '#5AA2E3', 'LineStyle', '--', 'LineWidth', 1.5);
% zero line
yline(0, ':k', 'LineWidth', 1.00)

xlim( [min(indepvar) max(indepvar)] ) 
grid("on")
box('on')

title('(b) EU rate', 'Interpreter', 'latex', 'Interpreter', 'latex')
legend('total', 'distribution', 'retention', 'EE reallocation', '', 'Interpreter','latex');
legend('boxoff')

if jj == 1
    ylim([ -1.75 0.75 ])
    yticks([-1.5 -1 -0.5 0 0.50])
    legend('location', 'southwest')
else
    legend('location', 'best')
end

hold off



% EE rate
subplot(2,2,3)
hold on

% total
plot(indepvar, ee.total, 'Color', '#0D16A0', 'LineWidth', 2.25);
% distribution
plot(indepvar, ee.skill + ee.job, 'Color', '#5AA2E3', 'LineStyle', '-.', 'LineWidth', 1.5);
% contact
plot(indepvar, ee.contact, 'Color', '#5AA2E3', 'LineStyle', ':', 'LineWidth', 1.5);
% selection
plot(indepvar, ee.selection, 'Color', '#5AA2E3', 'LineStyle', '--', 'LineWidth', 1.5);
% zero line
yline(0, ':k', 'LineWidth', 1.00)

xlim( [min(indepvar) max(indepvar)] ) 
grid("on")
box('on')

xlabel( x_label, 'Interpreter', 'latex' )
title('(c) EE rate', 'Interpreter', 'latex', 'Interpreter', 'latex')
legend('total', 'distribution', 'contact', 'selection', '', 'Interpreter','latex');
legend('boxoff')

if jj == 1
    ylim([ -0.80 0.15 ])
    yticks([-0.75 -0.5 -0.25 0])
    legend('location', 'southwest')
else
    legend('location', 'best')
end

hold off


% U rate
subplot(2,2,4)
hold on

% total
plot(indepvar, u.total, '-b', 'Color', '#0D16A0', 'LineWidth', 2.25);
% distribution
plot(indepvar, u.distr, 'Color', '#5AA2E3', 'LineStyle', '--', 'LineWidth', 1.5);
% transitions
plot(indepvar, u.trans, 'Color', '#5AA2E3', 'LineStyle', ':', 'LineWidth', 1.5);
% zero
yline(0, ':k', 'LineWidth', 1.00)

xlim( [min(indepvar) max(indepvar)] )
grid("on")

xlabel( x_label, 'Interpreter', 'latex' )
title('(d) Unemployment rate', 'Interpreter', 'latex', 'Interpreter', 'latex')
legend('total', 'distribution', 'conditional transition', '', 'Interpreter','latex');
legend('boxoff')
box('on')

if jj == 1
    ylim([ -1.5 1.5 ])
    legend('Location', 'northwest')
else
    legend('location', 'best')
end

hold off

% export figure
saveas(fig, [results_path,'/Figures/',filename], 'epsc')
saveas(fig, [results_path,'/Figures/',filename], 'png')

end


%%

disp('Decomposition: done. Check latex files.')
