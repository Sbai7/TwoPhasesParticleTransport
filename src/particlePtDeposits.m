function [particle,diff] = particlePtDeposits(Grid,V,particle,dt) 
% Solve the ode's of particle kinetics related to the elemental mechanism  
% of deposition in pore throats. 
%
% INPUTS:
% Grid              - Grid used for discretization 
% V                 - Structure of mid-cells edges velocities along x,y,and
%                     z directions
% particle          - instance of the particles subject to transport
%                     mechanisms
% dt                - Time interval over which to perform ode integration
%
% OUTPUTS:
% particle          - returns a modifield instance of the particles with  
%                     updated deposited concentrations in pore throats.
% diff              - cell-centered array of the net change in deposited 
%                     pore throats concentrations during this time step. It 
%                     is related to the change in mobile particle
%                     concentrations. 
%
% Author: M.A. Sbai, Ph.D.
%         BRGM (French Geological Survey) 
%         D3E  (Direction Eau, Environnement, Echotechnologies)
% 

Nx = Grid.Nx; Ny = Grid.Ny; Nz = Grid.Nz; 
N  = Grid.N;

% We need to interpolate the velocity accurately from cell-edges to get a
% mean value at the cell center => Just use Pollock's formula's for a
% regular/orthogonal cell to do this.

% Get the net velocity norm for the hydrodynamic rate expression
u = zeros(3,Nx,Ny,Nz);
u(1,1:Nx,:,:) = (V.x(1:Nx,:,:) + V.x(2:Nx+1,:,:))/2;
u(2,:,1:Ny,:) = (V.y(:,1:Ny,:) + V.y(:,2:Ny+1,:))/2;
u(3,:,:,1:Nz) = (V.z(:,:,1:Nz) + V.z(:,:,2:Nz+1))/2;
u_norm     = zeros(Grid.Nx,Grid.Ny,Grid.Nz);
u_norm(:,:,:) = sqrt( u(1,:,:,:).^2 + u(2,:,:,:).^2 + u(3,:,:,:).^2 );
u_norm = reshape(u_norm,N,1);

% Keep track of old deposited conc's.
C_pt_old = particle.C_pt;
%C_pt     = zeros(size(C_pt_old));

%--- Solve the kinetics ode.

% create an ode options structure
opts = odeset('AbsTol',1e-3,'RelTol',1e-2,'BDF','on','Stats','off');

% time span for integration 
tspan = [0,dt]; 

for i=1:N
   % initial value for the ode 
   y0 = C_pt_old(i);
   
   % shortcut to function handle to integarte with extra parameters
   fun = @(t,y) Rate_pt(t,y,particle.C(i),particle.alpha_pt,u_norm(i));
   
   % Sole the local ODE system with the stiff ode solver ODE15S with
   % backwards differentiation formula
   [~, buffer] = ode15s(fun, tspan, y0, opts);
   particle.C_pt(i,1) = buffer(size(buffer,1),1);
end

% Calculate rate 
diff = particle.C_pt - C_pt_old;

%--------------------------------------------------------------------------
% Nested rate law function for ODE integration 
%--------------------------------------------------------------------------
function ydot = Rate_pt(t, y, ...
                         C, ...
                         apt, ... 
                         u_norm)

ydot    = apt*u_norm.*C;           % Pore throat dep. rate
end

end