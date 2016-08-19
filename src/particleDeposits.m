function [C_dep,Rate]=particleDeposits(Grid,V,particle,C,C_dep,Cs,dt) 

% C_dep: input/output concentration of net deposited particles (could be 
%     negative) from colloidal and hydrodynamic release and deposition
%     kinetics.
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

net_u_norm = zeros(Grid.Nx,Grid.Ny,Grid.Nz);
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
C_dep_old = zeros(N,1); C_dep_old = C_dep;

% Solve the kinetics ode.
opts = odeset('abstol',1e-3,'reltol',1e-2,'stats','off');
for i=1:N
   [t, buffer] = daspk('Rate_d', [0 dt], C_dep_old(i), opts, ...
                  C(i), ... 
                  particle.alpha_h, particle.alpha_cl, particle.alpha_d, ... 
                  net_u_norm(i), net_conc(i), u_norm(i));
   C_dep(i) = buffer(size(buffer,1),1);
end

% Calculate rate 
Rate = zeros(N,1);
Rate = (C_dep - C_dep_old);
