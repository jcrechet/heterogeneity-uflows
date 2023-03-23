% Créchet J.
% created: december 2020; last updated March 2023.
% project: Heterogeneity in labor mobility and unemployment across
% countries
% Produce figure 1 and 2 in section 3 (baseline model)

%% Figure 1

load( "Results.mat" )

% rename tables for convenience
t1 = elasticities_UM;
t2 = elasticities_HM;
clearvars elasticities_uniform elasticities_hetero 

% parameter values to be plotted
llambda = [ t1.llambda, t2.llambda ];
sigma_x = t2.sigma_x;
bb = t1.bb;

% Figure
fig = figure;
hold on
fig.Position = [1000 250 1000 600]; 

% (a) UE elasticity
subplot('Position', [0.035 0.55 0.45 0.30])

plot( bb, t1.ue, "k:", bb, t2.ue_tightness, "b--", bb,  t2.ue_total, "b-", 'LineWidth', 1 )
xlim([0.55,0.95])
ylim([0 0.3])
title("(a) $| \, d \, \ln\Lambda_{UE} / d \, F \, |$ ", "Interpreter", "latex")
legend("UM elasticity", "HM elasticity, \textit{tightness} component", "HM total elasticity (\textit{tightness} + \textit{selection})", "Interpreter", "latex", Location="northwest")
grid("on")
xlabel("reservation wage ($b$)", "Interpreter", "latex")
yticks( [0 0.1 0.2 0.3] )

% (b) EU elasticity
subplot('Position', [0.535 0.55 0.45 0.30])

plot( bb, t1.eu, "k:", bb, t2.eu_retention, "b--", bb, t2.eu_total, "b-", 'LineWidth', 1 )
xlim([0.55,0.95])
ylim([0 0.3])
title("(b) $| \, d \, \ln\Lambda_{EU} / d \, F \, |$", "Interpreter", "latex")
legend("UM elasticity", "HM elasticity, \textit{retention} component", "HM total elasticity (\textit{retention} + \textit{distribution})", "Interpreter", "latex", Location="northwest")
grid("on")
xlabel("reservation wage ($b$)", "Interpreter", "latex")
yticks( [0 0.1 0.2 0.3] )

% (c) parameters
subplot('Position', [0.3 0.10 0.45 0.30])

plot( bb, llambda(:,1), "k:",  bb, llambda(:,2), "b--", bb, sigma_x.^2, "b-", 'LineWidth', 1 )
xlabel("reservation wage ($b$)", "Interpreter", "latex")

ylim([0 0.4])
yticks( [0 0.1 0.2 0.3 0.4] )

grid("on")
legend("$\lambda(b)$, UM model", "$\lambda(b)$, HM model", "$\sigma_x^2(b)$, HM model", "Interpreter", "latex", Location="northwest")
title("(c) Parameters", "Interpreter", "latex")

hold off

% save figures
saveas(fig, results_path+"\Figures\Figure_1", 'epsc')
saveas(fig, results_path+"\Figures\Figure_1", 'png')

%% Figure 2

% 'preferred' (illustrative) calibration
p.ddelta = 0.0039;
p.b = 0.7541;
p.llambda = 0.0837;
p.sigma_x = 0.4649;
m.eu = 0.0173;  

% compute benchmark (F=0)
p.F = 0;
[~, ~, ~, v0] = model_HM(p, true, []);

% compute counterfactual (F=3)
p.F = 3;
[~, ~, ~, v1]= model_HM(p, true, []);

% axis
mu_x = - p.sigma_x^2/2; % unconditional mean of match quality equal one
x_lb = exp(- 1 );       %exp( mu_x - 3.5*p.sigma_x );
x_ub = exp( mu_x + 3.5*p.sigma_x );
x = linspace(x_lb,x_ub,1000)';
log_x = log(x);

% unconditional distributions of match quality
gx = @( x ) lognpdf(x, mu_x, p.sigma_x);
Gx = @( x ) logncdf(x, mu_x, p.sigma_x);

% log 
k = @(x) log(x);

% variables
xR0 = v0.xR;
xR1 = v1.xR;
xhat0 = v0.xhat;
xhat1 = v1.xhat;
sx0 = v0.sx;
sx1 = v1.sx;
hx0 = v0.hx;
hx1 = v1.hx;


%% Figure 2, panel (a): Eql with no firing costs

tiledlayout(3,1)

% unconditional distribution
nexttile
plot( k(x), Gx(x), 'b-', 'LineWidth', 1 )

yline( Gx(xR0), 'b:',  'LineWidth', 1.5, ...
    'Label', '$G_x( \underline{x}_R )$', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )

xline( k(1), 'k-', 'LineWidth', 0.75 )


xline( k(xR0), 'b:',  'LineWidth', 1.5, ...
    'Label', '$\log( \underline{x}_R )$', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )

xline( k(xhat0), 'b:', 'LineWidth', 1.5, ...
    'Label', '$\log( \hat{x} )$', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )

axis([k(x_lb), k(x_ub), 0, 1.25])

ylabel('$G_x(x)$', 'Interpreter', 'latex'  )

% equilibrium distribution
nexttile

plot( k(x), hx0(x), 'b-', 'LineWidth', 1 )

xline( k(1), 'k-', 'LineWidth', 0.75 ) 
xline( k(xR0), 'b:',  'LineWidth', 1.5, ...
    'Label', '$ $', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )

xline( k(xhat0), 'b:', 'LineWidth', 1.5, ...
    'Label', '$ $', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )

axis([k(x_lb), k(x_ub),0 , 2.5])

ylabel('$h_x(x)$', 'Interpreter', 'latex'  )

% separation probability
nexttile

plot( k(x), log(sx0(x)), 'b-', 'LineWidth', 1 )

yline( log(m.eu), 'k:', 'LineWidth', 1,  'Label', '$\ln \big( \Lambda_{EU} \big)$', 'Interpreter', 'latex' )

xline( k(1), 'k-', 'LineWidth',  0.75)

xline( k(xR0), 'b:',  'LineWidth', 1.5, ...
    'Label', '$ $', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )

xline( k(xhat0), 'b:', 'LineWidth', 1.5, ...
    'Label', '$ $', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )
axis([k(x_lb),k(x_ub),-6,0])

xlabel('$\ln(x)$', 'Interpreter','latex')
ylabel('$\ln(s(x))$', 'Interpreter', 'latex'  )

% save figures
saveas(gcf, results_path+"\Figures\Figure_2a", 'epsc')
saveas(gcf, results_path+"\Figures\Figure_2a", 'png')


%% Figure 2, panel (b): Eql with high firing costs

tiledlayout(3,1)

% unconditional distribution
nexttile
plot( k(x), Gx(x), 'b-', 'LineWidth', 0.75 )

yline( Gx(xR0), 'b:',  'LineWidth', 1.5, ...
    'Label', '', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )

yline( Gx(xR1), 'r--',  'LineWidth', 1.5, ...
    'Label', '$G_x(x_R^\prime)$', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )

xline( k(1), 'k-', 'LineWidth', 0.75 )

xline( k(xR0), 'b:',  'LineWidth', 1.5, ...
    'Label', '$ $', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )

xline( k(xhat0), 'b:', 'LineWidth', 1.5, ...
    'Label', '$ $', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )

xline( k(xR1), 'r--',  'LineWidth', 1.5, ...
    'Label', '$ \log(\underline{x}_R^\prime)$', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left' )

xline( k(xhat1), 'r--', 'LineWidth', 1.5, ...
    'Label', '$ \log(\hat{x}^\prime)$', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left' )

axis([k(x_lb), k(x_ub), 0, 1.25])

ylabel(' ', 'Interpreter', 'latex'  )


% equilibrium distribution
nexttile

plot( k(x), hx0(x), 'b-', k(x), hx1(x), 'r-', 'LineWidth', 0.75 )

xline( k(1), 'k-', 'LineWidth', 0.75 ) 

xline( k(xR0), 'b:',  'LineWidth', 1.5, ...
    'Label', '$ $', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )
xline( k(xhat0), 'b:', 'LineWidth', 1.5, ...
    'Label', '$ $', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )

xline( k(xR1), 'r--',  'LineWidth', 1.5, ...
    'Label', '$ $', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )
xline( k(xhat1), 'r--', 'LineWidth', 1.5, ...
    'Label', '$ $', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )

axis([k(x_lb), k(x_ub), 0, 2.5])

ylabel(' ', 'Interpreter', 'latex'  )


% separation probability
nexttile

plot( k(x), log(sx0(x)), 'b-', k(x), log(sx1(x)), 'r-', 'LineWidth', 0.75 )

yline( log(m.eu), 'k:', 'LineWidth', 1,  'Label', '$\ln \big( \Lambda_{EU} \big)$', 'Interpreter', 'latex' )

xline( k(1), 'k-', 'LineWidth',  0.75)

xline( k(xR0), 'b:',  'LineWidth', 1.5, ...
    'Label', '$ $', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )
xline( k(xhat0), 'b:', 'LineWidth', 1.5, ...
    'Label', '$ $', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )

xline( k(xR1), 'r--',  'LineWidth', 1.5, ...
    'Label', '$ $', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )

xline( k(xhat1), 'r--', 'LineWidth', 1.5, ...
    'Label', '$ $', 'Interpreter', 'latex', 'LabelOrientation', 'horizontal' )

axis([k(x_lb), k(x_ub), -6, 0] )

xlabel('$\ln(x)$', 'Interpreter','latex')
ylabel(' ', 'Interpreter', 'latex'  )

% save figures
saveas(gcf, results_path+"\Figures\Figure_2b", 'epsc')
saveas(gcf, results_path+"\Figures\Figure_2b", 'png')