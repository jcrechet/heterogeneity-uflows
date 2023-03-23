function distance = calibration_secular(x, p, targets)


% initialize structure
p.calibration = true;
p.display_equil_iterations = false;
p.max_iter_eq = 50;

% param values
p.b = x(1);
p.chi = x(2);

% eval equil
[~, m] = equilibrium([], p, []);

% distance with targets
xx = [m.agg.b_w,   targets.b_w; ...
      m.agg.var_w, targets.var_w ];

% compute distance
distance = sum( abs( xx(:,1) - xx(:,2) ) ./ xx(:,2) );


end