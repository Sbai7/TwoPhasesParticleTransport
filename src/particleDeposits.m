function [particle,diff] = particleDeposits(Grid,V,particle,Cs,dt) 
% Solve the ode's of particle kinetics related to elemental mechanisms of 
% colloidal release from pore surfaces, hydrodynamic release from pore
% bodies, and deposition in pore surfaces. 
%
% INPUTS:
% Grid              - Grid used for discretization 
% V                 - Structure of mid-cells edges velocities along x,y,and
%                     z directions
% particle          - instance of the particles subject to transport
%                     mechanisms
% Cs                - Gridded array of cell-centered salinities used to
%                     trigger the colloidal release mechanism 
% dt                - Time interval over which to perform ode integration
%
% OUTPUTS:
% particle          - returns a modifield instance of the particles with  
%                     updated deposited concentrations in pore bodies.
% diff              - cell-centered array of the net change in deposited 
%                     pore bodies concentrations during this time step. It 
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

net_u_norm = u_norm;
for i=1:N
   net_u_norm(i) = max( net_u_norm(i) - particle.uc, 0 );
end

% Get net concentration for the colloidal rate expression
net_conc = zeros(N,1);
for i=1:N
   net_conc(i) = max( particle.Csc - Cs(i), 0);
end

% Keep track of old deposited conc's.
C_dep_old = particle.C_dep;

%--- Solve the kinetics ode.

% create an ode options structure
opts = odeset('AbsTol',1e-3,'RelTol',1e-2,'BDF','on','Stats','off');

% time span for integration 
tspan = [0,dt]; 

for i=1:N
   % initial value for the ode 
   y0 = C_dep_old(i);
   
   % shortcut to function handle to integarte with extra parameters
   fun = @(t,y) Rate_d(t,y,particle.C(i),particle.alpha_h,particle.alpha_cl, ...
       particle.alpha_d,net_u_norm(i), net_conc(i), u_norm(i));
   
   % Sole the local ODE system with the stiff ode solver ODE15S with
   % backwards differentiation formula
   [~, buffer] = ode15s(fun, tspan, y0, opts);
   particle.C_dep(i,1) = buffer(size(buffer,1),1);
end

% Calculate variation in deposited mass  
diff = particle.C_dep - C_dep_old;

%--------------------------------------------------------------------------
% Nested rate law function for ODE integration 
%--------------------------------------------------------------------------
function ydot = Rate_d(t, y, ...
                         C, ...
                         ah, acl, ad, ... 
                         net_u_norm, net_conc, u_norm)

ydot    = - ah*y.*net_u_norm ...    % Hydrodynamic rel. rate
          - acl*y.*net_conc  ...    % Colloidal rel. rate 
          + ad*u_norm.*C;           % Surface dep. rate
end

end