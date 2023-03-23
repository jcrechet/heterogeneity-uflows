function [llambda, elasticities] = model_UM(parameters, moments)

% Compute semi-elasticity of UE and EU wrt firing costs (assuming F=0)
% in the baseline model with uniform mobility (i.e., DMP model with
% ex-ante identical matches and idiosyncratic shocks---endogenous
% separations)

% arguments: parameters: structure with model parameter values; 
% moments: structure with targeted U, UE, and EU transition probabilities.

%%%
%%%

% rename structure arguments for convenience
p = parameters; 
m = moments; 
clear parameters moments


%%%
%%%


% 1. Solve for value of lambda (proabability of a productivity shock) 
% consistent with b and teh targeted EU rate (assuming uniform [0,1])

% specify separation threshold zr as a function of lambda given parameters and
% targeted eu rate
zr_lambda = @(ll) ( m.eu - p.ddelta ) ./ ( ( 1 - p.ddelta ) .* ll); 

% specify lambda as implicit function of other parameters (consistent with
% b and EU rate)
g = @(ll) zr_lambda(ll) - p.b + ( (p.bbeta*ll)/(1-p.bbeta*(1-ll))) * (1/2 - zr_lambda(ll)*( 1 - 1/2*zr_lambda(ll)) ); 

% fsolve to find root of function g
opts = optimoptions(@fsolve, 'Display', "off");
[llambda, fval, exit_flag] = fsolve( @(ll) g(ll), m.eu, opts);



%%%
%%%


% 3. if model can match data, compute elasticities
if ( exit_flag >= 1 && fval <= 1.00e-5 )

% initialize table
elasticities = table;

% 3.1. Compute surplus and z_r threshold given lambda

% compute zr associated with solution for lambda
zr = zr_lambda(llambda);

% deduce the value of the surplus
S0 = (p.z0-zr)/(1-p.bbeta*(1-llambda));


%%%
%%%


% 4. Compute semi-elasticities given lambda and z_r

% UE rate
d_ue = - ( (1-p.eeta) / p.eeta) * p.bbeta * m.ue * m.eu / ( 1 - p.bbeta * (1-m.eu) ) * S0^(-1);
elasticities.ue = abs( d_ue/m.ue );


% EU
d_eu = - ( 1 - p.bbeta) * ( 1 - p.bbeta * ( 1 - llambda ) ) / ( 1 - p.bbeta * ( 1 - m.eu)) * llambda;
elasticities.eu = abs( d_eu/m.eu ); 

else
    
    % if no solution, return nan
    llambda = nan;
    elasticities = table;
    elasticities.eu = nan;
    elasticities.ue = nan;
    
    disp(["no solution found !!, b = ", p.b])

end


%%%
%%%

end
