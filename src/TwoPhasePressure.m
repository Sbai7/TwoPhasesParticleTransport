function [P,V]=TwoPhasePressure(Grid,S,Fluid1,Fluid0,q,P,dt)

% When the last two arguments (P,dt) are supplied, calculations
% are assumed to be performed in transient conditions. Where 'P'
% is the solution at the previous time step and 'dt' is the current
% time step size.

% gravity = 9.81;

% Compute K*lambda(S) ... the viscous term
[M1,M0]=RelativePerm(S,Fluid1,Fluid0);
Mt = M1+M0;
KM = reshape([Mt,Mt,Mt]',3,Grid.Nx,Grid.Ny,Grid.Nz).*Grid.K;

% kz = reshape(Grid.K(3,:,:,:),Grid.Nx,Grid.Ny,Grid.Nz);
% Mg = zeros(Grid.Nx,Grid.Ny,Grid.Nz);
% Mg = reshape(Mw',Grid.Nx,Grid.Ny,Grid.Nz).*Fluid.rw + reshape(Mo',Grid.Nx,Grid.Ny,Grid.Nz).*Fluid.ro;
% KG = gravity.*Mg.*kz;

% Compute pressure and extract fluxes
if (nargin == 6)
   [P,V] = TpfaTwoPhase(Grid,KM,q,P,dt); % transient solver
else
   [P,V] = TpfaTwoPhase(Grid,KM,q);      % steady-state solver
% [P,V,rhs] = TpfaTwoPhase(Grid,KM,KG,q);
end