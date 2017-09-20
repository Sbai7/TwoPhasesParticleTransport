function [C,Flux]  = Diffusion(Grid,D,C,dt)
% Solve the three-dimensional diffusion equation with a tensor diffusion 
% coefficient.
%
% INPUTS:
% Grid              - Grid used for discretization 
% D                 - Diffusion tensor of the process 
% C                 - Concentration array at previous time step 
% dt                - Time interval over which to run the solver
%
% OUTPUTS:
% C                 - Array of new cell-centered concentrations
% Flux              - Diffusive mass flux at cellsd edges
%
% Author: M.A. Sbai, Ph.D.
%         BRGM (French Geological Survey) 
%         D3E  (Direction Eau, Environnement, Echotechnologies)
% 

% Compute transmissibilities by harmonic averaging.
Nx=Grid.Nx; Ny=Grid.Ny; Nz=Grid.Nz; N=Nx*Ny*Nz;
hx=Grid.hx; hy=Grid.hy; hz=Grid.hz;
L = D.^(-1);
tx = 2*hy*hz/hx; TX = zeros(Nx+1,Ny,Nz);
ty = 2*hx*hz/hy; TY = zeros(Nx,Ny+1,Nz);
tz = 2*hx*hy/hz; TZ = zeros(Nx,Ny,Nz+1);
TX(2:Nx,:,:) = tx./(L(1,1:Nx-1,:,:)+L(1,2:Nx ,:,:));
TY(:,2:Ny,:) = ty./(L (2,:,1:Ny-1,:)+L(2,:,2:Ny,:));
TZ(:,:,2:Nz) = tz./(L (3,:,:,1:Nz-1)+L(3,:,:,2:Nz));

% Assemble TPFA discretization matrix.
x1 = reshape(TX(1:Nx,:,:),N,1); x2 = reshape(TX(2:Nx+1,:,:),N,1);
y1 = reshape(TY(:,1:Ny,:),N,1); y2 = reshape(TY(:,2:Ny+1,:),N,1);
z1 = reshape(TZ(:,:,1:Nz),N,1); z2 = reshape(TZ(:,:,2:Nz+1),N,1);
DiagVecs = [-z2,-y2,-x2,x1+x2+y1+y2+z1+z2,-x1,-y1,-z1];
DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
A = spdiags(DiagVecs,DiagIndx,N,N);

% RHS vector
rhs = zeros(N,1);

% Add contribution from old solution 'C' to RHS vector and the system 
% matrix 
phiS = Grid.sat.*Grid.por./dt;
tterm = reshape(phiS,N,1);
A = A + spdiags(tterm,0,N,N);
C = reshape(C,N,1);
rhs = rhs + (tterm.*C);

% Solve linear system and extract diffusive cell-interface fluxes.
u = A\rhs;
C = reshape(u,Nx,Ny,Nz);
Flux.x = zeros(Nx+1,Ny,Nz);
Flux.y = zeros(Nx,Ny+1,Nz);
Flux.z = zeros(Nx,Ny,Nz+1);
Flux.x(2:Nx,:,:) = (C(1:Nx-1,:,:)-C(2:Nx,:,:)).*TX(2:Nx,:,:); % -D*dC/dx
Flux.y(:,2:Ny,:) = (C(:,1:Ny-1,:)-C(:,2:Ny,:)).*TY(:,2:Ny,:); % -D*dc/dy
Flux.z(:,:,2:Nz) = (C(:,:,1:Nz-1)-C(:,:,2:Nz)).*TZ(:,:,2:Nz); % -D*dC/dz 

C = reshape(C,N,1);

end
