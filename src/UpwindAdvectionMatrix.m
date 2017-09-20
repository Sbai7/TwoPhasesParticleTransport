function A = UpwindAdvectionMatrix(Grid,V,q)
% Generates the pure upwind advection matrix for a mass/heat/particle 
% transport process from a spatial distribution of the velocity field. 
%
% INPUTS:
% Grid              - Grid used for discretization 
% V                 - Structure of mid-cells edges velocities along x,y,and
%                     z directions
% q                 - Injection/Production flow rates 
%
% OUTPUT:
% A                 - Sparse matrix of the advection process
% 
% Author: M.A. Sbai, Ph.D.
%         BRGM (French Geological Survey) 
%         D3E  (Direction Eau, Environnement, Echotechnologies)
% 

Nx = Grid.Nx; Ny = Grid.Ny; Nz = Grid.Nz; N = Grid.N;

% extract production flow rates (i.e. negative part)
fp = min(q,0);

%--- We separate the velocity components along each direction into positive 
%    and negative parts

xn = min(V.x,0); 
yn = min(V.y,0); 
zn = min(V.z,0); 
xp = max(V.x,0); 
yp = max(V.y,0); 
zp = max(V.z,0); 

% reshape these arrays to column vectors 
xn = reshape(xn(1:Nx,:,:),N,1); 
yn = reshape(yn(:,1:Ny,:),N,1); 
zn = reshape(zn(:,:,1:Nz),N,1); 
xp = reshape(xp(2:Nx+1,:,:),N,1); 
yp = reshape(yp(:,2:Ny+1,:),N,1); 
zp = reshape(zp(:,:,2:Nz+1),N,1);

%--- The upwind advection sparse matrix stencil is composed from these  
%    seven diagonals. Production flow rates are add to the main diagonal  

DiagVecs = [zp,      yp, xp, fp+xn-xp+yn-yp+zn-zp, -xn, -yn,   -zn];
DiagIndx = [-Nx*Ny, -Nx, -1,         0,              1,  Nx, Nx*Ny];
A = spdiags(DiagVecs,DiagIndx,N,N);

end