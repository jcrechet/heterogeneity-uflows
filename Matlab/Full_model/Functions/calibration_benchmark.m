function [obj_rand_grid, x_rand_grid, obj0, x0] = calibration_benchmark(parameters,x_initial, run_rand_search, run_ga_minimization, run_ga_initial)

p = parameters;
clearvars parameters data

% Descrioption here

%% parameters for minimization

% random grid search options
nb_param = p.K;
pop_size = 100; % size of population and random grid

% simulated annealing minization options
if run_ga_minimization == 1

max_it = nb_param*10^3;
max_time = nb_param*(3600*24);
nb_elite = 5;
share_crossover = 0.5; 
nb_crossover = round( (pop_size - nb_elite)*share_crossover );
nb_mutation = pop_size-nb_elite-nb_crossover;
mut_var_scale = 5*10^(-2);  % mutation: eps \distras N(0,sigma), sigma = scale*(ub-lb);
pop_search_size = 10;       % size of population for pattern search
mod_patternsearch = 25;     % number of iterations before calling pattern search

end

%% objective function and bounds

% set appropriate options for estimation
p.calibration  = true;
p.equilibrium  = false;

% function to minimize
f = @(x) equilibrium(x, p, []);

% bounds
lb = p.lb;
ub = p.ub;

%% random grid search

if run_rand_search == 1

% draw matrix of random param vector
nb_draws = pop_size;

% initialize
x_rand_grid = zeros(nb_param,nb_draws);


% draw random values
for k = 1:nb_param    
    x_rand_grid(k,:) = lb(k) + (ub(k)-lb(k))*rand(1,nb_draws)';
end

% initial value (if provided)
if isempty(x_initial) == 0
x_rand_grid(:,1) = x_initial;
end


% evaluate model over random grid
% initial vector for storing value of objective
obj_rand_grid = zeros(nb_draws,1);

% loop over grid columns
parfor jj = 1:nb_draws
    obj_rand_grid(jj) = f(x_rand_grid(:,jj));
    disp(obj_rand_grid(jj))
end

% save intermediate results (random grid evaluation)
save('Results\calibration\outcomes_random_grid_search.mat', 'x_rand_grid', "obj_rand_grid")

%% Initialization for ga minimization + pattern search

% initialize grid

% load outcomes of random grid search
load('Calibration\outcomes_random_grid_search.mat', 'x_rand_grid', 'obj_rand_grid')

% initial population
[ obj_rand_grid_sorted, tmp_indexes ] = sort( obj_rand_grid ); % sort from best to worst
x_rand_grid_sorted = x_rand_grid(:,tmp_indexes);               % get associated parameters
x0 = x_rand_grid_sorted(:,1:pop_size);                         % use best points as initial guess in local search for min 
obj0 = obj_rand_grid_sorted(1:pop_size);

clear x_rand_grid obj_rand_grid x_rand_grid_sorted obj_rand_grid_sorted

% save initial grid for ga minimization
save('Calibration\ga_outcomes.mat', 'x0', 'obj0')

end


%% Run minimization algorithm

if run_ga_minimization == true

load('Calibration\ga_outcomes.mat', 'x0', 'obj0')


% option to evaluate grid on pre existing parameters
if run_ga_initial == false

disp( ' ')
disp('Evaluating objective function on a preexisting grid !')
disp( ' ')

% reset objective vector
obj0 = 0*obj0;

% parfor loop to evaluate function on pre-existing grid
parfor jj = 1:pop_size
   obj0(jj,1) = f(x0(:,jj)); 
end

disp('values: ')
disp( obj0 )

% update
save('Calibration\ga_outcomes.mat', 'x0', 'obj0')

end



% initialize loop
it = 0;
ga_timer = tic;

% initialize vector for new population outcomes in each iteration
x1 = x0;                                                       

while ( it < max_it && toc(ga_timer) < max_time )

%%%
    
% 1. elite
x1(:,1:nb_elite) = x0(:,1:nb_elite);

%%%

% 2. crossover

% form random couples of parents
parent = randperm(nb_crossover*2,nb_crossover*2)';
parent = reshape( parent, nb_crossover, 2 );

% crossover children
for i = 1:nb_crossover
    
    oomega = rand(nb_param,1);
    genes_1 = x0(:,parent(i,1));
    genes_2 = x0(:,parent(i,2));
    
    %x1(:,nb_elite+i) = (oomega<0.5).*genes_1 + (oomega>=0.5).*genes_2; 
    x1(:,nb_elite+i) = oomega.*genes_1 + (1-oomega).*genes_2; 

end

%%%

% 3. mutation
for i = 1:nb_mutation
    
    oomega = rand(nb_param,1);
    genes_1 = x0(:,parent(i,1));
    genes_2 = x0(:,parent(i,1));
    
    %x1(:,nb_elite+nb_crossover+i) = (oomega<0.5).*genes_1 + (oomega>=0.5).*genes_2 + randn(nb_param,1).*( mut_var_scale*(ub-lb) ) ;
    x1(:,nb_elite+nb_crossover+i) = oomega.*genes_1 + (1-oomega).*genes_2 + randn(nb_param,1).*( mut_var_scale*(ub-lb) ) ;
    x1(:,nb_elite+nb_crossover+i) = min( max( x1(:,nb_elite+nb_crossover+i), lb), ub);  
    
end

%%%

% 5. evaluate function over crossover and mutation populations
obj1 = zeros(pop_size,1);
obj1(1:nb_elite) = obj0(1:nb_elite);

parfor jj = nb_elite+1:pop_size
   obj1(jj,1) = f(x1(:,jj)); 
end

%%%

% 6. pattern search from a subset of the population
if (it>0 && mod(it,mod_patternsearch)==0)

   disp('run pattern search')

% descending sorting for selecting best points    
[ tmp_val, tmp_indexes ] = sort(obj1);
obj1 = tmp_val;
x1 = x0(:,tmp_indexes);

% run pattern search over current best population points
for jj = 1:pop_search_size
   [x1(:,jj), obj1(jj,1)] = patternsearch(f, x1(:,jj), [], [], [], [], lb, ub, [], p.pattern_options);
   % save
   save('Calibration\ga_outcomes.mat', 'x0', 'obj0', 'it', 'x1', 'obj1')
end
end

%%%

% 5. update
it = it+1;
[ obj0, tmp_indexes ] = sort(obj1); % fitness of new population
x0 = x1(:,tmp_indexes);             % new generation

%%%

% diagnostic
disp('generation')
disp(it)

disp('elite fitness')
disp(obj0(1:nb_elite))

disp('best genotype')
disp(x0(:,1))

% save workspace
save('Calibration\ga_outcomes.mat', 'x0', 'obj0', 'it')

end

if (it>max_it)
disp('max nb of iterations has been reached')
end

if (toc(ga_timer)>max_time)
disp('max time reached')
end

end



end