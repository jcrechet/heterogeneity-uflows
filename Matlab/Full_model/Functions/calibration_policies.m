function distance = calibration_policies(x, p, targets, step)


% initialize structure
p.display_equil_iterations = false;

% step 1: calibrate to policy targets
if step == 1

% param values
p.f = x(1);
p.b = x(2);
p.tax = x(3);

% eval equil
[~, m] = equilibrium([], p, []);

% distance with targets
xx = [m.agg.f_w ,  targets.f_w1; ...
      m.agg.b_w,   targets.b_w1; ...
      m.agg.tax_w, targets.tax_w1 ];


% step 1: calibrate matching efficiency to match UE rate difference
elseif step == 2

p.A = x(1);

% eval equil
[~, m] = equilibrium([], p, []);

xx = [m.agg.ue, targets.ue1];


% xx = [m.agg.ue, targets.ue1; ...
%      m.agg.eu, targets.eu1];

end

% compute distance
distance = sum( abs( xx(:,1) - xx(:,2) ) ./ xx(:,2) );


end