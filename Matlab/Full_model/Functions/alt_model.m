function distance = alt_model(x, p, targets, OJS)


% initialize structure
p.equilibrium = false;
p.calibration = true;


% values of parameters to recalibrate
p.A = x(1);
p.b = x(2);
p.lambda = x(3);

% case with OJS
if OJS == true
    p.s = x(4);
end


% evaluate equilibrium
[ ~, m ] = equilibrium([], p, []);


% compute distance
xx = [ m.agg.ue, targets.ue; m.agg.eu, targets.eu];


% case with OJS 
if OJS == true
    xx = [ xx; m.agg.ee, targets.ee  ];
end


% compute distance
distance = sum( abs( xx(:,1) - xx(:,2) ) ./ xx(:,2) );

if isnan( distance ) == true
    distance = inf;
end

end