function result_table = ee_variation(m0, m1, v0, v1, p)


% 1. Numerical computation of relative variations of EU rate and its components

% in: model structures m0 and m1 (moments) and v0 and v1 (equilibrium
% variables); index 0: benchmarck; 1: counterfactual, evaluated at the
% varied policy value; p: parameter structure
% out: a structure with a tables with results.

% Based on appendix E of paper.


%% 0. retrieve and construct required variables


% age range for aggregate outcomes (in model time units)
it1 = int64( ( p.age_range(1) - p.age_entry ) / p.periodlen );
it2 = ( p.age_range(2) - p.age_entry + 1 ) / p.periodlen - 1 ;


% (a) aggregate EU rate
EE0 = m0.agg.ee; 
EE1 = m1.agg.ee;


% (b) conditional EU rate
ee0 = v0.ee(:,it1:it2);
ee1 = v1.ee(:,it1:it2);


% (c) Joint distribution of skills and job type
h0 = v0.e(:,it1:it2) / sum( v0.e(:,it1:it2), "all" );  
h1 = v1.e(:,it1:it2) / sum( v1.e(:,it1:it2), "all" );


% (d) Employment marginal skill distribution

% Size of skill state-space dimension 
I1 = p.Ik;

% Size of job-type state-space
I2 = p.Ix*p.Iz;

% Size of age state space
T = length(h0(1,:));

% (i) benchmark 
tmp1 = reshape( h0, [I1, I2, T]);  % step 1: reshape matrix in 3 dim array form with dim 1: skill; dim 2: job; dim 3: age
tmp2 = sum( tmp1, 2);              % step 2: take sum over dim 2 (job) to get marginal skill pmf
tmp3 = permute( tmp2, [1, 3, 2]);  % step 3: permute so age becomes dim 2 
tmp4 = reshape( tmp3, I1, T );     % step 4: reshape 3 dim-array in matrix form
hw0 = repmat( tmp4, I2, 1 );       % step 5: adjust format to apply weights in upcoming calculations
clearvars tmp1 tmp2 tmp3 tmp4

% (ii) counterfactual
tmp1 = reshape( h1, [I1, I2, T]);
tmp2 = sum( tmp1, 2);
tmp3 = permute( tmp2, [1, 3, 2]);
tmp4 = reshape( tmp3, I1, T );
hw1 = repmat( tmp4, I2, 1 );
clearvars tmp1 tmp2 tmp3 tmp4

% (c) Job distribution conditional on skill
hj0 = h0 ./ hw0; hj0(h0 == 0) = 0;
hj1 = h1 ./ hw1; hj1(h1 == 0) = 0;

% (d) Contact
p0 = v0.ptheta * v0.s_e(it1:it2);
p1 = v1.ptheta * v1.s_e(it1:it2);

% (e) selection
G0 = ee0 ./ p0;
G1 = ee1 ./ p1;

% cleanup
clearvars m0 m1 v0 v1 p


%% 1. Compute actual change
actual = ( EE1 - EE0 ) / EE0 ;


%% some checking
% EE0_ = nansum( hj0.*hw0.*p0.*G0, "all" );
% EE1_ = nansum( hj1.*hw1.*p1.*G1, "all" );
% actual_ =  ( EE1_ - EE0_ ) / EE0_ ;



%% 2. Compute weights

% 2.1. all components
pphi = ( (h1 + h0).*( ee0 + ee1 ) ./ (4*EE0) );

% 2.2. distribution subcomponents
ggamma1 = ( hw0 + hw1 ).*( hj0 + hj1 ) ./ ( 2*(h0 + h1) );

% 2.3. conditional-transition subcomponents
ggamma2 = ( p0 + p1 ).*( G0 + G1 ) ./ ( 2*(ee0 + ee1) );



%% 3. Compute components

% 2.1. skills
tmp = ( pphi.*ggamma1 ).*( 2*( hw1 - hw0 ) ./ ( hw1 + hw0 ) );
tmp( pphi == 0 ) = 0;
skill = sum( tmp, "all" );

% 2.2. job
tmp = ( pphi.*ggamma1 ).*( 2*( hj1 - hj0 ) ./ ( hj1 + hj0 ) );
tmp( pphi == 0 ) = 0;
job = sum( tmp, "all" );

% 2.3. contact
tmp = ( pphi.*ggamma2 ).*( 2*( p1 - p0 ) ./ ( p1 + p0 ) );
tmp( pphi == 0 ) = 0;
contact = sum( tmp, "all" );

% 2.4. selection
tmp =  ( pphi.*ggamma2 ).*( 2*( G1 - G0 ) ./ ( G1 + G0 ) );
tmp( pphi == 0 ) = 0;
selection = sum( tmp, "all" );

% 2.5. total
total = skill + job + contact + selection;
residual = actual - total;



%% 4. Contruct table
result_table = table(actual, total, residual, skill, job, contact, selection);


end