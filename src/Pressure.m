function [P,V] = Pressure(Grid,Fluid,q,P,dt)
% Drives the single phase pressure solver in steady-state and transient
% modes. 
%
% INPUTS:
% Grid              - Grid used for discretization 
% Fluid             - Fluid properties 
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

%--- Compute k/mu array
k_mu = zeros(3,Nx,Ny,Nz);

if isscalar(Fluid.viscosity)    % uniform fluid viscosity
    
    k_mu = Grid.K / Fluid.viscosity;
    
else                            % array of cell-centered fluid viscosity
    Fluid.viscosity = Fluid.viscosity(:);
    assert(length(Fluid.viscosity)==N);
    for i=1:3
        k_mu(i,:,:,:) = Grid.K(i,:,:,:) ./ ...
                        reshape(Fluid.viscosity,1,Nx,Ny,Nz);
    end
    
end

%--- Compute pressure and velocities
if (nargin > 4)
    
    % call two-point flux approximation solver in transient mode
   [P,V] = TPFA(Grid,k_mu,q,P,dt);
   
else
    
	% call two-point flux approximation solver in steady-state mode
    [P,V] = TPFA(Grid,k_mu,q);
    
end

end