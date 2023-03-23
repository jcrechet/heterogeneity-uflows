function [distance, moments, elasticities, variables] = model_HM(parameters, compute_elasticities, x )


% Compute semi-elasticity of UE and EU wrt firing costs (assuming F=0)
% in the baseline model with heterogeneous mobility (i.e., DMP model with
% stochastric matching and idiosyncratic shocks)

% arguments: parameters: structure with model parameter values; 
% moments: structure with targeted U, UE, and EU transition probabilities.
% distance: distance with targeted data


%%%
%%%


% rename structure function arguments
p = parameters; 
clear parameters moments


%%%
%%%


% 1. parameters and functional forms/distributions
if isempty(x) == false
    p.llambda = x(1);
    p.sigma_x = x(2);
end

% unconditional mean of match quality equal one
mu_x = - p.sigma_x^2/2;

% specify distribution
gx = @(x) lognpdf(x,mu_x,p.sigma_x); % pdf of match quality
Gx = @(x) logncdf(x,mu_x,p.sigma_x); % cdf 

%%%
%%%


% 2. separation probability as a function of x given parameters

% 2.1. solve for separation z reservation threshold given x

% root coefficients for quadratic equation (solve for z such that surplus equal zero)
a = @(x) p.bbeta*p.llambda*x/(2*(1-p.bbeta));
b = @(x) x;
c = @(x) p.bbeta*p.llambda*(x-2*p.b)/(2*(1-p.bbeta)) - p.b + (1-p.bbeta*(1-p.llambda))*p.F ;

% compute zR as a function of x
zR_x = @(x) (c(x) >= 0).*(0) ...
   + (a(x)+b(x)+c(x) <= 0).*(1) ...
   + (c(x) < 0 & a(x)+b(x)+c(x) > 0).*( - b(x)./(2*a(x)) + ((b(x).^2-4*(a(x).*c(x))).^(1/2)) ./ (2*a(x)) );

% 2.2. solve for hiring threshold
if p.F==0
    xR = p.b;
else
    tmp_f = @(x) (1-zR_x(x))*x - (1-p.bbeta*(1-p.llambda))*p.F;
    x0 = p.b+p.bbeta*p.llambda;
    xR = fsolve(tmp_f, x0);
end

% 2.3. deduce separation probability as a function of match quality
sx = @(x) (x >= xR) .* ( p.ddelta + (1-p.ddelta) * p.llambda * zR_x(x) );


%%%
%%%


% 3. Compute model moments

% Match quality steady-state equilibrium distribution
gtilde = @(x)  gx(x)./sx(x);
tmp = integral( gtilde, xR, inf);
hx = @(x) (x>=xR).*gtilde(x)/tmp;

% (i) aggregate EU rate
eu = ( 1 - Gx(xR) ) / tmp;

% (ii) separation rate by tenure 
h_tn = cell(120,1);     % cell for match quality distribution by tenure
s_tn = zeros(120,1);    % vector for separation profile

% initial period (i.e., first period)
h_tn{1} = @(x) (gx(x)./(1-Gx(xR))).*(x>=xR); 
integrand = @(x) sx(x).*h_tn{1}(x);
s_tn(1) = integral(integrand,xR,inf);

% subsequent periods
for itn = 2:120
    h_tn{itn} = @(x) ((1-sx(x))./(1-s_tn(itn-1))).*h_tn{itn-1}(x);
    integrand = @(x) sx(x).*h_tn{itn}(x);
    s_tn(itn) = integral(integrand,xR,inf);
end

% compute yearly separation profile
m = (1:120)';
y = floor( m./12 )';
sr_y = zeros(10,1);
for iy = 0:9
    sr_y(iy+1) = mean(s_tn(y==iy));
end

% store in structure
moments = struct;
moments.eu = eu;
moments.sr_y = sr_y;


% 4. Model distance to targets

% model ratio eu rate zero year of tenure to eu rate bw 5 to 9 year  
eu_ratio = sr_y(1)/( mean( sr_y(6:10) ) );

% distance
tmp1 = [ eu,  0.0173; eu_ratio, 0.0308/0.0061 ] ;
tmp2 = abs( tmp1(:,1) - tmp1(:,2) ) ./ tmp1(:,2);
distance = mean( tmp2 );

clearvars tmp1 tmp2

moments.eu_ratio = eu_ratio;

%%%
%%%


% 4. Elasticities
if compute_elasticities == true

    %%%

    % 4.1. UE rate

    % inaction cutoff
    xhat = ( 1 - p.bbeta * ( 1 - p.llambda )) / ( p.bbeta * p.llambda)* 2*( p.b - (1-p.bbeta)*p.F ) ;

    % hiring surplus function
    S0 = @(x) ( x >= xR & x < xhat ) .* ( (1 - zR_x(x)).*x ./ (1-p.bbeta*(1-p.llambda)) ) ...
        + ( x >= xhat ) .* (  (1-p.bbeta*(1 - p.llambda))^(-1)*( 1 + p.bbeta*p.llambda / (2*(1-p.bbeta)) )*x - p.b/(1-p.bbeta) + p.F  ) ;

    % expected hiring surplus (expected profits, \gamma=0)
    ES0 = integral( @(x) S0(x).*gx(x), xR, inf);

    % derivative of the increase in expected separation costs
    % wrt to F
    fun1 = @(x) zR_x(x) ./  ( 1 - p.bbeta*( 1 - p.llambda*zR_x(x) ) ) .* gx(x);
    expected_cost = p.bbeta*p.llambda * integral( fun1, xR, inf );

    % Elasticity of tightness component
    ue_tightness = (1-p.eeta)/p.eeta * expected_cost * ES0^(-1);

    % Elasticity of selection component
    ue_selection = p.bbeta*p.llambda * gx(xR)/(1-Gx(xR));

    % Total
    ue_total = ue_tightness + ue_selection;


    %%%

    % 4.2. EU rate

    % integrand for the 'direct', retention effect
    fun1 = @(x) (x>=xR & x<=xhat).*(1-p.bbeta*(1-p.llambda))./(1-p.bbeta*(1-p.llambda*zR_x(x))).*(x.*sx(x)).^(-1).*hx(x);

    % integrand for the 'composition' effect
    fun2 = @(x) (1-zR_x(x)).*hx(x);

    % Elasticity of retention component
    eu_retention = (1-p.bbeta)*p.llambda*integral( fun1, xR, xhat);

    % Elasticity of composition component
    eu_composition = p.bbeta*p.llambda*gx(xR)/(1-Gx(xR))*integral( fun2, xR, inf);

    % total
    eu_total = eu_retention + eu_composition;

    %%%

    % stack in a table
    elasticities = table(ue_tightness, ue_selection, ue_total, eu_retention, eu_composition, eu_total);

    %%%

    % model variables
    variables = struct;
    variables.gx = gx;
    variables.hx = hx;
    variables.xR = xR;
    variables.xhat = xhat;
    variables.sx = sx;

%%%
%%%

else

    elasticities = [];
    variables = [];

end

end