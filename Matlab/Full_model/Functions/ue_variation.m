function result_table = ue_variation(m0, m1, v0, v1, p)

% 1. Numerical computation of relative variations of the UE rate and its components

% in: model structures m0 and m1 (moments) and v0 and v1 (equilibrium
% variables); index 0: benchmarck; 1: counterfactual, evaluated at the
% varied policy value; p: parameter structure
% out: a structure with tables with the results.

% Based on appendix E of the paper.


%% 0. retrieve required variables

% age range for aggregate outcomes (in model time units)
it1 = int64( ( p.age_range(1) - p.age_entry ) / p.periodlen );
it2 = ( p.age_range(2) - p.age_entry + 1 ) / p.periodlen - 1 ;

% aggregate UE rate
UE0 = m0.agg.ue; 
UE1 = m1.agg.ue;

% conditional UE rate
ue0 = v0.ue(:,it1:it2);
ue1 = v1.ue(:,it1:it2);

% unemployment skill distribution
lu0 = v0.u(:,it1:it2)/sum( v0.u(:,it1:it2), "all" );
lu1 = v1.u(:,it1:it2)/sum( v1.u(:,it1:it2), "all" );

% tightness
theta0 = v0.theta;
theta1 = v1.theta;

% proba of contact per time unit
ptheta0 = v0.ptheta;
ptheta1 = v1.ptheta;

% Elasticity of matching
Eta = p.Eta;

% Efficiency of matching
A0 = ptheta0 / ( theta0^(1-Eta) );
A1 = ptheta1 / ( theta1^(1-Eta) );

% search
s0 = v0.s_u(:,it1:it2);
s1 = v1.s_u(:,it1:it2);

% selection
G0 = ue0 ./ ( ptheta0*s0 );
G1 = ue1 ./ ( ptheta1*s1 );

% cleanup
clearvars m0 m1 v0 v1 p


%% 1. Compute actual change
actual = ( UE1 - UE0 ) / UE0 ;

% implied change (some checking)
% UE0_ = sum( (s0*ptheta0).*G0.*lu0, "all");
% UE1_ = sum( (s1*ptheta1).*G1.*lu1, "all");
% actual_ = ( UE1_ - UE0_ ) / UE0_ ;


%% 2. Compute weights

% all components
pphi = ( lu0 + lu1 ).*( ue0 + ue1 ) ./ (4*UE0) ;
pphi( isnan(pphi) ) = 0;

% all components in transition rate
pphi2 =  (1/4) * ( ( theta0^(1-Eta) + theta1^(1-Eta) )*(A0*s0+A1*s1).*(G0+G1) ./ ( ue0 + ue1 ) );
pphi2( isnan(pphi2) ) = 0;

% tightness
alpha1 = pphi2.* ( (2/3) * ( ( (A0*s0).*G0 + (A1*s1).*G1 ) ./ ( ( (A0*s0)+(A1*s1) ).*( G0+G1 ) ) + 1) );
alpha1( isnan( alpha1 ) ) = 0;

% search
alpha2 = pphi2.* ( (2/3) * ( ( theta0^(1-Eta) * G0 + theta1^(1-Eta) * G1 ) ./ ( ( theta1^(1-Eta) + theta0^(1-Eta) ).*( G0 + G1 ) ) + 1) );
alpha2( isnan( alpha2 ) ) = 0;

% selection
alpha3 = pphi2.* ( (2/3) * ( ( ptheta0.*s0 + ptheta1.*s1 ) ./ ( ( ptheta0+ptheta1 ).*( s0+s1 ) ) + 1) );
alpha3( isnan( alpha3 ) ) = 0;


%% 2. Compute components

% 2.1. Skill 

% compute elements
tmp = pphi.*( lu1 - lu0 )./( ( lu1 + lu0 )/2 );

% apply zero weigths where required
tmp(pphi==0) = 0;

% take the sum
skill = sum( tmp, "all" );

% 2.2. tightness
tmp = pphi.*alpha1.*( theta1^(1-Eta) - theta0^(1-Eta) )./( ( theta1^(1-Eta) + theta0^(1-Eta) )/2 );
tmp(pphi==0) = 0;
tightness = sum( tmp, "all" );


% 2.3. search
tmp = pphi.*alpha2.*( A1*s1 - A0*s0 )./( ( A1*s1 + A0*s0 )/2 );
tmp(pphi==0) = 0;
search = sum( tmp, "all" );

% 2.4. selection
tmp = pphi.*alpha3.*( G1 - G0 )./( ( G1 + G0 )/2 );
tmp(pphi==0) = 0;
selection = sum( tmp, "all" );

% total computed
total = skill + tightness + search + selection;
residual = actual - total; 


%% 3. make table
result_table = table(actual, total, residual, skill, tightness, search, selection);


end