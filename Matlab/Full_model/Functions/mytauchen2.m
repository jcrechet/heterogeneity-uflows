% myTauchen2: adapted from Jan Hannes Lang code "mytauchen" (written
% 6.7.2010)
% available at: https://sites.google.com/site/janhanneslang/programs
% (consulted 15.11.2022)---description is below
% Project: Heteroegeneity in labor mobility and unemployment flows across
% countries, Créchet (2022)
% this function generates a common discretized skill grid for (1) the
% unemployment state (2) the unemployment state
% the two states being different in terms of the mean for the associated
% AR(1) processes.

function [s, Pi0, Pi1] = mytauchen2(N,mu0,mu1,rho,sig,m)

s       = zeros(N,1);
s(1)    = mu0/(1-rho) - m*sqrt(sig^2/(1-rho^2));
s(N)    = mu1/(1-rho) + m*sqrt(sig^2/(1-rho^2));
step    = (s(N)-s(1))/(N-1);

cdf_normal = @(x) 0.5 * erfc(-x/sqrt(2));


for i=2:(N-1)
   s(i) = s(i-1) + step; 
end

% transition matrix for U state
Pi0      = zeros(N,N);
for j = 1:N
    for k = 1:N
        if k == 1
            Pi0(j,k) = cdf_normal((s(1) - mu0 - rho*s(j) + step/2) / sig);
        elseif k == N
            Pi0(j,k) = 1 - cdf_normal((s(N) - mu0 - rho*s(j) - step/2) / sig);
        else
            Pi0(j,k) = cdf_normal((s(k) - mu0 - rho*s(j) + step/2) / sig) - ...
                      cdf_normal((s(k) - mu0 - rho*s(j) - step/2) / sig);
        end
    end
end

% transition matrix for E state
Pi1      = zeros(N,N);
for j = 1:N
    for k = 1:N
        if k == 1
            Pi1(j,k) = cdf_normal((s(1) - mu1 - rho*s(j) + step/2) / sig);
        elseif k == N
            Pi1(j,k) = 1 - cdf_normal((s(N) - mu1 - rho*s(j) - step/2) / sig);
        else
            Pi1(j,k) = cdf_normal((s(k) - mu1 - rho*s(j) + step/2) / sig) - ...
                      cdf_normal((s(k) - mu1 - rho*s(j) - step/2) / sig);
        end
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: [s, Pi] = mytauchen(mu,rho,sig,N)
%
% This function discretizes a continuous AR(1) process by using the method
% proposed by Tauchen (1986). The AR(1) process takes the following form:
% y(t) = mu + rho*y(t-1) + eps(t), where eps ~ N(0,sig^2)
% Parts of the code are taken from the function tauchen.m written by Martin
% Flodén.
%
% INPUTS
%   mu:     scalar, intercept term of the AR(1) process
%   rho:    scalar, AR-coefficient
%   sig:    scalar, standard deviation of innovations
%   N:      scalar, number of grid points for the discretized process
%
% OUTPUTS
%   s:      column vector of size Nx1, contains all possible states in ascending order
%   Pi:     matrix of size NxN, contains the transition proabilities. Rows
%           are current state and columns future state
%   
% Author:   Jan Hannes Lang
% Date:     6.7.2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%