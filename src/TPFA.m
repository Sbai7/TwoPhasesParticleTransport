function [P,V] = TPFA(Grid,K,q,P,dt)
% Solves 3D pressure equation under steady-state or transient conditions. 
%
% INPUTS:
% Grid              - Grid used for discretization 
% K                 - Apparent permeability tensor (i.e. permeability is 
%                     multiplied by another parameter such as cell centered
%                     mobility of a given phase)
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
hx = Grid.hx; hy = Grid.hy; hz = Grid.hz;

%--- Compute transmissibilities by harmonic averaging
L  = K.^(-1);
tx = 2*hy*hz/hx; 
ty = 2*hx*hz/hy; 
tz = 2*hx*hy/hz; 
TX = zeros(Nx+1,Ny,Nz); TY = zeros(Nx,Ny+1,Nz); TZ = zeros(Nx,Ny,Nz+1);
TX(2:Nx,:,:) = tx./(L(1,1:Nx-1,:,:)+L(1,2:Nx,:,:));
TY(:,2:Ny,:) = ty./(L(2,:,1:Ny-1,:)+L(2,:,2:Ny,:));
TZ(:,:,2:Nz) = tz./(L(3,:,:,1:Nz-1)+L(3,:,:,2:Nz));

%--- Assemble TPFA discretization matrix
x1 = reshape(TX(1:Nx,:,:),N,1); x2 = reshape(TX(2:Nx+1,:,:),N,1);
y1 = reshape(TY(:,1:Ny,:),N,1); y2 = reshape(TY(:,2:Ny+1,:),N,1);
z1 = reshape(TZ(:,:,1:Nz),N,1); z2 = reshape(TZ(:,:,2:Nz+1),N,1);
DiagVecs = [-z2,-y2,-x2,x1+x2+y1+y2+z1+z2,-x1,-y1,-z1];
DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
A        = spdiags(DiagVecs,DiagIndx,N,N);

% In transient mode (i.e. when P and dt arguments are supplied) add
% contribution from old solution 'P' to RHS vector and the system matrix
if (nargin == 5) 
   phic  = Grid.compr.*Grid.por./dt;
   tterm = reshape(phic,N,1);
   A     = A + spdiags(tterm,0,N,N);
   P     = reshape(P,N,1);
   q     = q + (tterm.*P);
end

%--- Solve resulting sparse linear system of equations 
P = A\q;
P = reshape(P,Nx,Ny,Nz);

%--- Calculate mid-cells edges velocities using Darcy's law
V.x = zeros(Nx+1,Ny,Nz); V.y = zeros(Nx,Ny+1,Nz); V.z = zeros(Nx,Ny,Nz+1);
V.x(2:Nx,:,:) = (P(1:Nx-1,:,:)-P(2:Nx,:,:)).*TX(2:Nx,:,:);
V.y(:,2:Ny,:) = (P(:,1:Ny-1,:)-P(:,2:Ny,:)).*TY(:,2:Ny,:);
V.z(:,:,2:Nz) = (P(:,:,1:Nz-1)-P(:,:,2:Nz)).*TZ(:,:,2:Nz);

end
