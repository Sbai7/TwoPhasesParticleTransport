function [C_pt,Rate]=particlePtDeposits(Grid,V,particle,C,C_pt,dt) 

% C_pt: input/output concentration of particles deposited in pore throats.
% dt: time step over which ode's integration is performed.


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
C_pt_old = zeros(N,1); C_pt_old = C_pt;

% Solve the kinetics ode.
opts = odeset('abstol',1e-3,'reltol',1e-2,'stats','off');
for i=1:N
   [t, buffer] = daspk('Rate_pt', [0 dt], C_pt_old(i), opts, ...
                  C(i), ... 
                  particle.alpha_pt, ... 
                  u_norm(i));
   C_pt(i) = buffer(size(buffer,1),1);
end

% Calculate rate 
Rate = zeros(N,1);
Rate = (C_pt - C_pt_old);
