function result_table = eu_variation(m0, m1, v0, v1, p)


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
EU0 = m0.agg.eu; 
EU1 = m1.agg.eu;


% (b) conditional EU rate
eu0 = v0.eu(:,it1:it2);
eu1 = v1.eu(:,it1:it2);


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


% (d) Probability of potential separation shock

% probality of exogenous separation
Delta  = p.Delta;

% transiton matrix for employment state
Pi_e = v0.Pi_e;

% preallocate arrays
Gz0 = zeros( I1*I2, T ); 
Gz1 = zeros( I1*I2, T ); 

% policy decision indicator: separation conditional on realized period state and no contact
tmp0 = ( v0.S(:,it1+1:it2+1) <= 0 );     
tmp1 = ( v1.S(:,it1+1:it2+1) <= 0 );   

% loop over age to compute probability of separation conditional on
% predetermined period state
for it = 1:T
    Gz0(:,it) = Delta + (1-Delta)*Pi_e*tmp0(:,it);
    Gz1(:,it) = Delta + (1-Delta)*Pi_e*tmp1(:,it);
end

clearvars tmp0 tmp1 


% (e) Probability of separation conditional on potential separation
% shock
eu0_nr = eu0./Gz0;
eu1_nr = eu1./Gz1;


% cleanup
clearvars m0 m1 v0 v1 p


%% 1. Compute actual change
actual = ( EU1 - EU0 ) / EU0 ;


%% some checking
% EU0_ = nansum( hj0.*hw0.*eu0_nr.*Gz0, "all" );
% EU1_ = nansum( hj1.*hw1.*eu1_nr.*Gz1, "all" );
% actual_ =  ( EU1_ - EU0_ ) / EU0_ ;



%% 2. Compute weights

% 2.1. all components
pphi = ( (h1 + h0).*( eu0 + eu1 ) ./ (4*EU0) );

% 2.2. distribution subcomponents
ggamma1 = ( hw0 + hw1 ).*( hj0 + hj1 ) ./ ( 2*(h0 + h1 ) );

% 2.3. conditional-transition subcomponents
ggamma2 = ( Gz0 + Gz1 ).*( eu0_nr + eu1_nr ) ./ ( 2*(eu0 + eu1 ) );



%% 3. Compute components

% 2.1. skills
tmp = ( pphi.*ggamma1 ).*( 2*( hw1 - hw0 ) ./ ( hw1 + hw0 ) );
tmp( pphi == 0 ) = 0;
skill = sum( tmp, "all" );

% 2.2. job
tmp = ( pphi.*ggamma1 ).*( 2*( hj1 - hj0 ) ./ ( hj1 + hj0 ) );
tmp( pphi == 0 ) = 0;
job = sum( tmp, "all" );

% 2.3. retention
tmp = ( pphi.*ggamma2 ).*( 2*( Gz1 - Gz0 ) ./ ( Gz1 + Gz0 ) );
tmp( pphi == 0 ) = 0;
retention = sum( tmp, "all" );

% 2.4. reallocation
tmp =  ( pphi.*ggamma2 ).*( 2*( eu1_nr - eu0_nr ) ./ ( eu1_nr + eu0_nr ) );
tmp( pphi == 0 ) = 0;
reallocation = sum( tmp, "all" );

% 2.5. total
total = skill + job + retention + reallocation;
residual = actual - total;



%% 4. Contruct table
result_table = table(actual, total, residual, skill, job, retention, reallocation);


end