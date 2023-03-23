function [distance, moments, variable_struct, timer_struct] = equilibrium(x, parameters, ptheta_initial)


% Function name: equilibrium.m
% Project: Heterogeneity in labor mobility and unemployment flows across countries
% Créchet, 2022
% created 2020 - last updated January 2022

% Compute a numerical equilibrium solution

% IN:

% x: vector of internally calibrated parameter values
% parameters: structure with parameter settings
% for solution algorithm, simulations, model, etc. (note that the values in the vector
% x prevail on those specified in parameters)
% ptheta_initial: scalar with initial value for p(theta)

% OUT:

% distance: measure of a distance (a mmetric value) between targeted empirical and computed moments.
% moments: computed moments
% variable_struct: structure containing model's equilibrium variables 
% (value functions, LM tightness, conditional transition rates, etc)
% timer_struct: a structure with timer results for different blocks of code (VFI,
% equilibrium distribution, etc).



%% Structure renaming and initialization

% rename parameter input struct
p = parameters; clearvars parameters

% initialize structure to store equilibrium variables
vb = struct;

% structure for timer
timer_struct = struct;

% start timer
timer_start = tic;


%% Assign internally calibrated parameter values from input vector x

if isempty(x) == false

% non-work utility
p.b     = x(p.ib);

% matching efficiency
p.A    = x(p.iA);

% probability of idiosyncratic shock to match output
p.lambda = x(p.ilambda);

% efficiency of employment search relative to unemployment
p.s  = x(p.is);              

% AR(1) variance of log ability (disturbance)
p.sigma_a = x(p.isigma_a);   

% AR(1) persistence of log ability
p.rho_a = x(p.irho_a);

% ability growth rates  
p.mu1 = x(p.imu1);

% variance of log permanent match quality
p.sigma_x = x(p.isigma_x);     

% search cost linear term
p.c0    = x(p.ic0); 

% passive effort threshold
p.slb   = x(p.islb); 

end

% normalization: match quality has mean 1
p.mu_x = -p.sigma_x^2/2; 


%% 1. Specify search and output functions

% (i) Search-cost function

% benchmark case: endogenous search intensity
if p.searchInt == 1 
    
    % cost function
    cost = @(s) ( p.c0^(1/p.c1) / (1+p.c1) )*max(s-p.slb, 0).^(1+p.c1);       
    
    % optimal search effort as a function of marginal return to search
    sopt = @(mr) min(  p.slb + max(mr,0).^(1/p.c1) / p.c0, 1 );      

% condition for counterctual: shutting down endogenous search intensity
else                
    
    % cost is zero
    cost   = @(s) 0;
    
    % search is exog., constant and equal one
    sopt   = @(s) 1;

end


% (ii) Match-output function
f = @(k,z,x) z .* exp( p.chi * ( log(x) + log(k) ) );



%% 2. Discretize state space and construct grids

tic

% (i) skills

% benchmark: heterogenous skills
if p.Ik > 1

    % use function "mytauchen2" to create TM and vector of skills
    m = 2;
    [ log_k, Pk0, Pk1 ] = mytauchen2(p.Ik, p.mu0, p.mu1, p.rho_a, p.sigma_a, m);
    k_vec = exp( log_k );
    
    % specify the indice for k-state of newborn worker
    Elog_a0 = 0; % mean of log equal 0
    [~, ia0] = min( abs( Elog_a0 - log_k ) );

    % cleanup
    clearvars log_k Elog_a0 m 

% counterfactual: shut down skill het.
else
    
    Pk0 = 1;
    Pk1 = 1;
    k_vec = 1;
    ia0 = 1;

end


% (ii) permanent match quality

if p.Ix > 1

    % use function "tauchen" to get the pmf vector for distribution match quality at hiring (with
    % rho = 0)
    m = 2.5; rho = 0;
    [log_x, Ptmp] = tauchen(p.Ix, p.mu_x, rho, p.sigma_x, m);
    
    % pmf vector is equal to any one of the columns of the TM (constructed using rho = 0) 
    gx = Ptmp(1,:)';  

    % adjust gx so it sum up to one
    gx = gx + (1-sum(gx))/p.Ix;

    % match-quality employment continuation TM: identity (match quality is constant within matches)
    Px = eye(p.Ix);  
    x_vec = exp(log_x); %

    % cleanup
    clearvars Ptmp log_x

else

    x_vec = 1;
    gx = 1;
    Px = 1;

end


% (iii) stochastic match output component

% vector for discretized uniform [0,1]
z_vec = linspace(0,1,p.Iz)';

% TM: iid shocks with probability = lambda
Pz = (1-p.lambda)*eye(p.Iz) + p.lambda*(1/p.Iz);


% (iv) Output grid and full state-space vector

% state-space grid size
Iu = p.Ik;             % unemployment state size
Ie = p.Ik*p.Iz*p.Ix;   % employment state size

% create "ND grids" using vectors for discretized state variables
% (grids in N=3 dimensional space)
[k_ndgrid, z_ndgrid, x_ndgrid] = ndgrid(k_vec, z_vec, x_vec);

% vectorize the N=3 dim grids to compute VF in vector form 
k_grid = reshape(k_ndgrid, Ie, 1);
z_grid = reshape(z_ndgrid, Ie, 1);
x_grid = reshape(x_ndgrid, Ie, 1);

% deduce the ND match output grid in vector form
y_grid = f(k_grid, z_grid, x_grid);

% indices for z-state in new potential job
% assuming z_0 = 1 (Mortensen, Pissarides; 1994)
iz0 = p.Iz;                  

% identify the indexes associated with z=z_0 
% on the full ND vector grids 
iz0_grid = ( z_grid == z_vec(iz0) );


% (v) Exogenous state transition matrices for VFI and equilibrium distributions

% Employment 
Pi_e = kron( kron(Px, Pz), Pk1);

% sparsing to save memory and fasten computation
Pi_e  = sparse( Pi_e );

% Unemployment
Pi_u  = sparse( Pk0 );

% timer
timer_struct.grid = toc;


% (vi) additional matrices for computing steady-state distribution

% transition matrix for worker state conditional on eu transition
Pi_eu = sparse( kron( ones(p.Iz*p.Ix,1), Pk1) );

% transition matrices for new match sampling conditional on realized
% state

% new z sampling pmf
gz = zeros(p.Iz,1);
gz(iz0) = 1;

% unemployment sampling matrix conditional on realized human capital state
Gu = sparse( kron( kron( gx, gz )', eye( p.Ik ) ) );

% employment sampling matrix conditional on realized human capital state
Ge = sparse( kron( ones(p.Iz*p.Ix,1), Gu ) );

% matrix B for uu and eu human capital transitions (for taking sum of probabilities cond. on human capital) 
B = zeros(Ie, Iu);
for i = 1:Ie
for j = 1:Iu
    B(i,j) = ( k_grid(i) == k_vec(j) ) ;
end
end
B = sparse(B);

clearvars k_ndgrid z_ndgrid x_ndgrid k_grid z_grid x_grid

%% 3. Value function backward resolution and equilibrium tightness

% start timer for equilibrium computation
equilibrium_timer = tic;

% (i) Initialization

% Initial guess for tightness 
if isempty(ptheta_initial) == 1
    theta_ = 1;
    ptheta_ = p.A;
else
    ptheta_ = ptheta_initial;
    theta_ = ( ptheta_ / p.A )^( 1 / (1-p.Eta) );
end

% Initialize distance measure for equilibrium fixed point computation (theta)
dist_eq = inf; 

% initialize number of iteration tracker for equil. computation
iter_eq = 0;   


% (ii) Value function terminal values

% time horizon
T = p.T;

% Unemployment
U(1:Iu, T) = p.b;                     

% Surplus
S(1:Ie, T) = ( (1-p.tax)*y_grid - p.b + p.f );

% Wage
W(1:Ie, T) = p.Gamma*(1-p.tax)*y_grid + p.Gamma*p.f + (1-p.Gamma)*p.b;

% Optimal search
s_u(1:Iu, T) = p.slb;
s_e(1:Ie, T) = p.slb;


% Start loop for iteration over equilbrium tigthness values
while dist_eq > p.crit_eq && iter_eq <= p.max_iter_eq

    % case with partial equilibrium (calibration): exit loop after first
    % iteration
    if p.equilibrium == 0
        dist_eq = 0;    
    end

    % update iteration tracker
    iter_eq = iter_eq + 1;

    % update value of theta and ptheta
    theta = theta_;
    ptheta = ptheta_;
    
    % timer for vf computation
    tic

    % (iii) Compute value functions backward given ptheta

    for it = T-1:(-1):1
            
            % (a) Unemployment

            % Expectation conditional on a contact and t+1 current worker
            % realized state
            Snew = reshape( max( S(iz0_grid, it+1) - p.f, 0), Iu, p.Ix); % array for potential new match
            EV = p.Gamma * Snew * gx;

            % Optimal search intensity
            MR = p.Beta * ptheta * EV;         % marginal return of search effort in unemployment
            s_u(:, it) = sopt(MR);             % optimal search effort

            % Compute the value of unemployment
            U(:, it) = p.b - cost( s_u(:, it) ) ...
                     + p.Beta*( ptheta*s_u(:, it).*( Pi_u * EV ) + Pi_u * U(:, it+1) );

            % Reservation wage when in employment
            wr = repmat( U(:, it) - p.Beta * Pk1 * U(:, it+1), p.Iz*p.Ix, 1 );


            % (b) Surplus function

            % Continuation surplus
            
            % cont. surplus in case of no contact and t+1 realized current match state
            Sstay = max( S(:, it+1), 0);

            % cont. surplus change conditional on contact and t+1 realized current match state
            Snew1 = repmat( Snew, [p.Iz*p.Ix, 1] ) ;
            tmpS  = max( Snew1 - Sstay, 0 );

            % Probability of quit conditional on a contact and t+1 realized current match state
            Proba_quit = ( tmpS > 0 ) * gx;

            % Worker marginal gain conditional on contact t+1 realized current match state
            Delta_W = p.Gamma * ( tmpS * gx );

            % Employer marginal loss cond. on contact and realized
            % match state
            Delta_J = - (1-p.Gamma) * ( Proba_quit .* Sstay );

            % Next-period surplus expectated change cond. on contact
            ES1 = Pi_e * ( Delta_W + Delta_J );

            % Next-period surplus conditional on no contact
            ES0 = Pi_e * Sstay;

            % Optimal search intensity
            MR = p.Beta * (1 - p.Delta) * p.s*ptheta * Pi_e * ( Delta_W + Delta_J + Proba_quit * p.f );
            s_e(:, it) = sopt(MR);

            % probability of a contact with an outside firm
            p_e = p.s*ptheta*s_e(:, it);  

            % compute surplus
            S(:, it) =  (1-p.tax)*y_grid - cost( s_e(:, it) )       ...                           % current period total return
                        + p.Beta*(1-p.Delta)*( p_e.*ES1 + ES0 )     ...                           % next-period expected discounted value
                        - wr + ( 1-p.Beta*(1-p.Delta)*( 1 - p_e.* ( Pi_e*Proba_quit ) ) ) * p.f;  % outside option 

            % wage
            W(:,it) =  p.Gamma*(1-p.tax)*y_grid + (1-p.Gamma)*cost( s_e(:,it) ) + (1-p.Gamma)*wr ...
                      + p.Gamma*( 1-p.Beta*(1-p.Delta)*( 1 - p_e.* ( Pi_e*Proba_quit ) ) ) * p.f ...
                      + p.Beta*p.Gamma*(1-p.Gamma)*p_e.* Pi_e *( Delta_J + Delta_W );

    end
    
    % cleanup
    clearvars Snew EV MR wr Snew1 tmpS Proba_quit Delta_W Delta_J ES0

    timer_struct.vf = toc;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % (iv) steady-state distribution 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % start timer for ss distribution
    tic

    % (a) policy decision TM
    
    % unemployment hiring
    Gamma_u = ( S >= p.f ); 
    
    % employment continuation
    Gamma_e = ( S >= 0 );

    % matrices for employment continuation conditional on a contact
    
    % initialize: TM conditional on realized match state and a contact
    if p.s > 0
    Gamma_ee = cell(T,1);

    % initialize: TM conditional on realized match state, a contact joint with EE
    % reallocation
    Gamma_ee_hat = cell(T,1);
    end
    
    % loop over age
    f = p.f;
    for it = 1:T
        
        % joint prob of employment continuation and EE reallocation cond.
        % on contact
        tmp1 = sparse( Ge .* ( S(:,it)' - f > max( S(:,it), 0 ) ) );
        
        % joint prob of staying in employment and no EE reallocation cond.
        % on contact
        tmp2 = sparse( ( 1 - sum( tmp1, 2 ) ) .* ( S(:,it) >= 0 )  );
        
        if p.s > 0

        % Gamma_hat
        Gamma_ee_hat{it} = tmp1 ;
        
        % Gamma
        Gamma_ee{it} = tmp1 + diag( tmp2 );

        end

    end
    clearvars tmp1 tmp2

    % probability of contact vector
    p_u = ptheta*s_u;
    p_e = ptheta*p.s*s_e;

    % (b) initialization

    % initialize vectors
    u = zeros(Iu, T); 
    u(ia0, 1) = 1/T;
    e = zeros(Ie, T);

    % vector for transition probabilities by state
    ue = zeros(Iu,T);
    eu = zeros(Ie,T);
    ee = zeros(Ie,T);

    % initialize cell for TM across employment state
    Puu = cell(T, 1);
    Pue = cell(T, 1);
    Peu = cell(T, 1);
    Pee = cell(T, 1);

    % cell for TM conditional on EE move (for EE rate and tenure
    % distribution)
    Pee_hat = cell(T, 1);
    

    % (c) loop over age to compute age-spe TM
    for it = 1:T-1
    
        % 1. employment to employment

        % probability of a contact
        p_cont = p_e(:,it);
        
        if p.s > 0
            Pee{it} = (1-p.Delta) ...
                      * ( (1-p_cont) .* ( Pi_e .* Gamma_e(:,it+1)' ) + p_cont .* ( Pi_e * Gamma_ee{it+1} ) );
        else
            Pee{it} = (1-p.Delta) * ( Pi_e .* Gamma_e(:,it+1)' );
        end
        % NOTE: Pi_e .* Gamma_e(:,it+1)' = Pi_e * diag( Gamma_e(:,it+1) );

        % 2. employment to unemployment
        Peu{it} = Pi_eu - Pee{it}*B;

        % 3. unemployment to employment
        Pue{it} = p_u(:,it) .* Pi_u * ( Gu .* Gamma_u(:,it+1)' );
        % Note: Pi_u .* Gamma_u(:,it+1)' = Pi_e * diag( Gamma_e(:,it+1) );
        
        % 4. unemployment to unemployment
        Puu{it} = Pi_u - Pue{it}*B;

        % 5. employer to employer
        % compute auxiliary matrix Pee_hat using Gamma_ee_hat
        if p.s > 0
            Pee_hat{it} = ( (1-p.Delta) * p_cont .* ( Pi_e * Gamma_ee_hat{it+1} ) );
        else
            Pee_hat{it} = zeros(Ie, Ie);
        end

        % sparse( (1-p.Delta) .* p_e(:,it) .* Pee_c_move );

        % 6. conditional transitions
        ue(:,it) = sum( Pue{it}, 2 );
        eu(:,it) = sum( Peu{it}, 2 ); 
        ee(:,it) = sum( Pee_hat{it}, 2 ); % job to job

    end

    clearvars p_cont 
    
    % (d) compute distribution backward
    for it = 2:T
        u(:,it) = u(:,it-1)'*Puu{it-1} + e(:,it-1)'*Peu{it-1} ;
        e(:,it) = u(:,it-1)'*Pue{it-1} + e(:,it-1)'*Pee{it-1} ;
    end

    timer_struct.ss_distribution = toc;


    % (e) profits of a vacancy
    
    % initialize
    EJ0 = 0; 

    % loop over experience
    for it = 1:T-1

        % mass of job seekers given ability
        u_t = s_u(:,it) .* u(:,it);
        e_t = p.s*s_e(:,it) .* e(:,it);

        % probability of contact with an unemployed worker
        prob_u = sum(u_t) / (sum(e_t) + sum(u_t));

        % sampling probability conditional on meeting with an unemployed worker
        if prob_u > 0
            h_u_t =  s_u(:,it) .* u(:,it) / sum(u_t);
        else
            h_u_t = 0*u_t;
        end

        % sampling probability conditional on meeting with an employed worker
        if (1-prob_u) > 0
            h_e_t =  p.s * s_e(:,it) .* e(:,it) / sum(e_t);
        else
            h_e_t = 0*e_t;
        end

        % employer profits, meeting unemployed worker
        Snew = reshape( max( S(iz0_grid, it+1) - p.f, 0), Iu, p.Ix); % array for potential new match
        EJ0_u = (1-p.Gamma) * ( (Pk0*h_u_t)' * (Snew * gx) );
        
        % employer profits, meeting with employed worker
        
        % transition matrices for state space
        Scurrent = max( S(:, it+1), 0);
        Snew1    = repmat( Snew, [p.Iz*p.Ix, 1] );

        EJ0_e = (1-p.Gamma) * ( (Pi_e*h_e_t)' * ( ( (Snew1 > Scurrent) .* Snew1 ) * gx ) );

        % update unconditional profits
        EJ0 = EJ0 + ( prob_u * EJ0_u + (1-prob_u) * EJ0_e ) * (1/(T-1));

    end

    
    % (f) compute tightness/vacancy posting cost if partial equilibrium
    if p.equilibrium == 0

    vb.c_v = p.A^(1/(1-p.Eta)) * ptheta^(-p.Eta/(1-p.Eta)) * p.Beta * EJ0;

    % (g bis) compute and update tightness and ptheta
    else
        
        % update theta: mid-point bw guess and value implied by updated
        % distribution and value functions
        theta_ = (1/2)*( ( p.A * p.Beta * EJ0 / p.c_v )^( 1/p.Eta ) + theta );
        
        % update ptheta to compute updated surplus and distributions
        ptheta_ = min( p.A * theta_^(1-p.Eta), 1);
        
        % measure distance between updated and initial ptheta
        dist_eq = abs( ptheta_ - ptheta );

        % display fixed-point search progress..
        if ( p.display_equil_iterations == true )
        disp('Iteration: ')
        disp(iter_eq)
        disp('theta:')
        disp(theta)
        disp('ptheta:')
        disp(ptheta)
        disp('ptheta, new iteration')
        disp(ptheta_)
        disp('distance')
        disp(dist_eq)
        end
        
        % check whether convergence occured
        if dist_eq < p.crit_eq 

            % diagnostic
            if p.display_equil_final == true
            disp('An equilibrium has been found!')
            end

            % store the computed fixed-point ptheta
            vb.ptheta = ptheta_;

            % equilibrium distributions
            vb.u = u;
            vb.e = e;

            % surplus
            vb.S = S;

            % search effort
            vb.s_u = s_u;
            vb.s_e = s_e;

            % tightness
            vb.theta = theta;

            % expected profits
            vb.EJ0 = EJ0;

            % transitions
            vb.eu = eu;
            vb.ue = ue;
            vb.ee = ee;

            % store transition matrix for employment state (useful for some
            % numerical decompositions)
            vb.Pi_e = Pi_e;

        end

    end

end

% Case with no convergence of tightness
if ( p.equilibrium == true && iter_eq >= p.max_iter_eq )
    disp('No equilibrium found for given tolerance; distance:')
    disp(dist_eq)
end

% checks transition matrices are well defined
% check that all entries are nonnegative and that matrices sum up to one
check_noneg   = zeros(T-1,1);
check_sumup_u = zeros(Iu,T-1);
check_sumup_e = zeros(Ie,T-1);

for it = 1:T-1
    check_noneg(it)     = 1 - (sum(Puu{it}<-eps,"all") + sum(Pue{it}<-eps,"all") + sum(Peu{it}<-eps,"all") + sum(Pee{it}<-eps,"all") );
    check_sumup_u(:,it) = sum(Puu{it},2) + sum(Pue{it},2);
    check_sumup_e(:,it) = sum(Peu{it},2) + sum(Pee{it},2);
end

% matrix diagnose
if sum( check_noneg < 1, "all" ) > 0
disp('transition matrices has negative elements')
end

if sum( abs(1 - check_sumup_u), "all" ) + sum( abs(1 - check_sumup_e), "all" ) > 10-6
    disp('matrix might not sum up to one')
end

clearvars check_noneg check_sumup_u check_sumup_e
clearvars U S s_u s_e Snew Snew1 

timer_struct.equilibrium = toc(equilibrium_timer);
variable_struct = vb;


%% 4. Moments for calibration


% initialize structure for storing model computed moments
m = struct;

% (a) aggregate statistics (age 20 to 60 in model time units) 
it1 = int64( ( p.age_range(1) - p.age_entry ) / p.periodlen );
it2 = ( p.age_range(2) - p.age_entry + 1 ) / p.periodlen - 1 ;

% create table  
m.agg = table;

% unemployment distribution
lu = u(:,it1:it2) / sum( u(:,it1:it2), "all" );

% employment distribution
he = e(:,it1:it2) / sum( e(:,it1:it2), "all" );

% aggregate flows
m.agg.ue = sum( ue(:,it1:it2) .* lu, "all" ); 
m.agg.eu = sum( eu(:,it1:it2) .* he, "all" ); 
m.agg.ee = sum( ee(:,it1:it2) .* he, "all" ); 

% unemployment rate
m.agg.urate = m.agg.eu / ( m.agg.eu + m.agg.ue ); 

% mean wage
m.agg.mean_W = sum( W(:,it1:it2) .* he, [1, 2] );

% replacement ratio 
m.agg.b_w = p.b / m.agg.mean_W;

% log wage mean
mean_w = sum( log(W(:,it1:it2)) .* he, [1, 2] );
m.agg.mean_w = mean_w;

% log wage variance
m.agg.var_w = sum( ((log(W(:,it1:it2)) - mean_w).^2) .* he, [1, 2] );
clearvars mean_w

% total output
y_tmp = repmat(y_grid, [1, T]); 
m.agg.Y = sum( y_tmp.*e, [1, 2] );

% match productivity
m.agg.mean_Y = m.agg.Y / sum( e, "all" );

% firing costs
m.agg.f_w = p.f / m.agg.mean_W;


% tax expressed in proportion of average (gross) wage

% average tax collected on the match
mean_tax = p.tax * m.agg.mean_Y;

% average pre-tax wage
mean_w_pre_tax = m.agg.mean_W + mean_tax;

% implied tax in proportion to pre-tax wage
m.agg.tax_w = mean_tax / mean_w_pre_tax;


% (b) experience profiles: UE/EU transitions and cumulative wage growth

% table
m.age = table;

% vector for experience in years
exp_year = floor( (1:T)*p.periodlen );

% preallocate age vectors
L = T*p.periodlen;
m.age.eu = zeros(L,1);
m.age.ue = zeros(L,1);
m.age.w = zeros(L,1);

% loop over age
for iy = 0:L-1

    % unemploynent by age
    u_tmp = u(:,exp_year==iy); e_tmp = e(:,exp_year==iy);
    u_yr = sum( u_tmp, [1, 2] ); e_yr = sum( e_tmp, [1, 2] );
    m.age.u(iy+1) = u_yr / ( e_yr + u_yr );

    % distributions cond. on age
    u_tmp = u_tmp./u_yr;
    e_tmp = e_tmp./e_yr;

    % transition rates
    ue_tmp = ue(:,exp_year==iy);
    eu_tmp = eu(:,exp_year==iy);

    % wage
    w_tmp = log( W(:,exp_year==iy) );

    % profiles
    m.age.ue(iy+1) = sum( ue_tmp.*u_tmp, [1, 2] );
    m.age.eu(iy+1) = sum( eu_tmp.*e_tmp, [1, 2] );
    m.age.w(iy+1)  = sum( w_tmp.*e_tmp, [1, 2] );

end

clearvars u_tmp e_tmp u_yr e_yr

% normalize log wage wrt wage after 1 year of experience
m.age.w = m.age.w - m.age.w(2);

% age
m.age.age = ( p.age_entry:p.age_entry+T*p.periodlen-1 )';

% age-specific transitions
% retrieve tmp vectors
u_ = m.age.u;
e_ = 1 - u_;
ue_ = m.age.ue;
eu_ = m.age.eu;
age_ = m.age.age;

% age 20-24
ind_tmp = ( age_ >= 20 & age_ <= 24 );
m.agg.ue_2024 = ( u_(ind_tmp) / sum(u_((ind_tmp))) )'*ue_(ind_tmp);
m.agg.eu_2024 = ( e_(ind_tmp) / sum(e_((ind_tmp))) )'*eu_(ind_tmp);

% age 25-29
ind_tmp = ( age_ >= 25 & age_ <= 29 );
m.agg.ue_2529 = ( u_(ind_tmp) / sum(u_((ind_tmp))) )'*ue_(ind_tmp);
m.agg.eu_2529 = ( e_(ind_tmp) / sum(e_((ind_tmp))) )'*eu_(ind_tmp);

% ue and eu rate, 30-39
ind_tmp = ( age_ >= 30 & age_ <= 39 );
m.agg.ue_3039 = ( u_(ind_tmp) / sum(u_((ind_tmp))) )'*ue_(ind_tmp);
m.agg.eu_3039 = ( e_(ind_tmp) / sum(e_((ind_tmp))) )'*eu_(ind_tmp);

% ue and eu rate, 40-49
ind_tmp = ( age_ >= 40 & age_ <= 49 );
m.agg.ue_4049 = ( u_(ind_tmp) / sum(u_((ind_tmp))) )'*ue_(ind_tmp);
m.agg.eu_4049 = ( e_(ind_tmp) / sum(e_((ind_tmp))) )'*eu_(ind_tmp);

% age 50-54
ind_tmp = ( age_ >= 50 & age_ <= 54 );
m.agg.ue_5054 = ( u_(ind_tmp) / sum(u_((ind_tmp))) )'*ue_(ind_tmp);
m.agg.eu_5054 = ( e_(ind_tmp) / sum(e_((ind_tmp))) )'*eu_(ind_tmp);

clearvars u_ e_ ue_ eu_ ee_ age_

%% 4b. Compute model distance with targeted statistics

if p.calibration == 1

    % vector moments
    xx = [ m.agg.ue, 0.3048; m.agg.eu, 0.0171; m.agg.ee, 0.0229; m.agg.var_w, .3180; ...
           m.agg.ue_2024, 0.3344; m.agg.ue_2529, 0.3195; m.agg.ue_3039, 0.3120; m.agg.ue_4049, 0.2952; m.agg.ue_5054, 0.2697; ...
           m.agg.eu_2024, 0.0340; m.agg.eu_2529, 0.0214; m.agg.eu_3039, 0.0168; m.agg.eu_4049, 0.0138; m.agg.eu_5054, 0.0120 ]; 
    
    % vector of absolute difference
    yy = abs((xx(:,1)-xx(:,2)))./xx(:,2);

    % take the mean
    distance = mean( yy ) ;

else

    distance = [];

end

%% 5. Compute additional distributions

if p.calibration == 0

    % (a) Compute unemployment distribution
    
    % initialize matrix for u distribution
    L = 25;                                   % up to two years (24 months)
    u_dur = zeros(Iu, T, L*(12*p.periodlen)); % dim 1: state, dim 2: age; dim 3: uduration

    % compute initial value: age zero (newborn)
    it = 1;
    u_dur(ia0, it, 1) = 1/T;

    % initial value: age 1, 2,...
    for it = 2:T

        u_dur(:, it, 1) = Peu{it-1}'*e(:,it-1);
    
    end

    % compute subsequent values

    % loop over age
    for it = 2:T

        % loop over u duration
        for id = 2:L*(12*p.periodlen)

            u_dur(:, it, id) = Puu{it-1}'*u_dur(:, it-1, id-1);
        
        end
        
    end


    % (b) Compute tenure distribution

    % initialize array
    L = 21;                             % tenure max: up to 21 years
    e_tn = zeros(Ie, T, L/p.periodlen); % dim 1: state; dim 2: age; dim 3: tenure

    % preallocate cell for TM joint with probability of e continuation
    Pee_ten = cell(T,1);
    
    % loop over age to compute age-specific matrix
    for it = 1:T

        % tenure ee stay matrix
        Pee_ten{it} = Pee{it}-Pee_hat{it};
    
    end

    for it = 2:T

        % initial mass
        e_tn(:, it, 1) = Pue{it-1}'*u(:,it-1) + Pee_hat{it-1}'*e(:,it-1);

        % subsequent
        for ie = 2:L/p.periodlen
            e_tn(:, it, ie) = Pee_ten{it-1}'*e_tn(:,it-1, ie-1);
        end

    end
    

    %% Compute additional moments 

    % (a) Age profiles: EE transitions and log-wage variance

    % preallocate vectors
    L = T*p.periodlen;
    m.age.ee = zeros(L,1);
    m.age.var_w = zeros(L,1);

    % loop over age
    for iy = 0:L-1

        % employment mass by age
        e_tmp = e(:,exp_year==iy);
        e_yr = sum( e_tmp, [1, 2] );

        % distribution cond. on age
        e_tmp = e_tmp./e_yr;

        % transition vector
        ee_tmp = ee(:,exp_year==iy);
        
        % wage-variance vector
        wmean = m.age.w(iy+1);
        w_tmp = ( log( W(:,exp_year==iy) ) - wmean ).^2;
        
        % profiles
        m.age.ee(iy+1) = sum( ee_tmp.*e_tmp, [1, 2] );
        m.age.var_w(iy+1)  = sum( w_tmp.*e_tmp, [1, 2] );

    end

    clearvars e_tmp e_yr

    % normalize wrt wage with one year of experience
    m.age.var_w = m.age.var_w - m.age.var_w(2);


    % (b) Unemployment duration UE profiles
    
    % table 
    m.udur = table;

    % preallocate
    L = 24;
    m.udur.ue = zeros(L,1);

    % vectorize ue
    ue_tmp = reshape( ue, [], 1);
    
    % vector for conditionning on month
    dur_month = (0:L); % to be adjusted if change in time unit

    for im = 0:L-1
        
        u_tmp = sum( u_dur(:, :, im == dur_month ), 3);
        u_tmp = reshape( u_tmp, [], 1);
        u_tmp = u_tmp / sum( u_tmp );
        m.udur.ue( im + 1 ) = u_tmp'*ue_tmp;

    end


    % (c) Tenure EU and EE profiles
    
    % table
    m.ten = table;

    % initialize 
    L = 21;
    m.ten.eu = zeros(L,1);
    m.ten.ee = zeros(L,1);

    % vector for conditionning on year
    ten_year = floor( (1:L/p.periodlen)*p.periodlen );

    % vector for separation conditional on state
    ee_tmp = reshape( ee, [], 1);
    eu_tmp = reshape( eu, [], 1);

    for iy = 0:L-1

    % vector for distribution
    e_tmp = sum( e_tn(:, :, iy == ten_year ), 3);
    e_tmp = reshape( e_tmp, [], 1 );
    e_tmp = e_tmp/sum(e_tmp);

    % compute transitions conditional on tenure
    m.ten.eu(iy+1) = e_tmp'*eu_tmp;
    m.ten.ee(iy+1) = e_tmp'*ee_tmp;

    end

end

clearvars Puu Pue Pee Peu

% rename moments structure
moments = m;

% timer
timer_struct.total = toc(timer_start);

end



