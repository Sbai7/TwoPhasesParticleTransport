function [P,V] = TwoPhasePressure(Grid,S_nw,Fluid1,Fluid0,q,P,dt)
% Drives the two-phase pressure solver in steady-state and transient modes. 
%
% INPUTS:
% Grid              - Grid used for discretization 
% S_nw              - Array of cell-centered fluid1 (non-wetting) phase 
%                     saturation 
% Fluid1            - Fluid properties of the non-wetting (injected) phase 
% Fluid0            - Fluid properties of the wetting (residant) phase
% q                 - Injection/Production flow rates 
% P  (optional)     - Pressure computed at previous time step 
% dt (optional)     - Time interval over which to run the transient solver
%
% When the last two arguments (P,dt) are supplied, calculations are assumed
% to be performed in transient conditions.
%
% OUTPUTS:
% P                 - Array of cell-centered pressure
% V                 - Structure of mid-cells edges velocities along x,y,and
%                     z directions
%
% Author: M.A. Sbai, Ph.D.
%         BRGM (French Geological Survey) 
%         D3E  (Direction Eau, Environnement, Echotechnologies)
% 

Nx = Grid.Nx; Ny = Grid.Ny; Nz = Grid.Nz; N = Grid.N;

%--- Compute the product K*M(S) where M is the total mobility 

% calculate mobilities of fluid1 and fluid0. M1 and M0, respectively.
[M1,M0] = RelativePerm(S_nw,Fluid1,Fluid0);

% total mobility
M_total = M1+M0;

% compute the product K*M
KM = reshape([M_total,M_total,M_total]',3,Nx,Ny,Nz).*Grid.K;

%--- Compute pressure and velocities
if (nargin > 6)
    
    % call two-point flux approximation solver in transient mode
    [P,V] = TPFA(Grid,KM,q,P,dt);

else
   
    % call two-point flux approximation solver in steady-state mode
    [P,V] = TPFA(Grid,KM,q);
   
end

end