function [P,V]=Pres(Grid,Fluid,q,P,dt)

% When the last two arguments (P,dt) are supplied, calculations 
% are assumed to be performed in transient conditions. Where 'P'  
% is the solution at the previous time step and 'dt' is the current
% time step size.

% Compute k/mu array 
k_mu = (Grid.K)./Fluid.viscosity;

% Compute pressure and extract fluxes
if (nargin == 5)
   [P,V] = TPFA(Grid,k_mu,q,P,dt); % transient solver
else 
   [P,V] = TPFA(Grid,k_mu,q);      % steady-state solver
end   
