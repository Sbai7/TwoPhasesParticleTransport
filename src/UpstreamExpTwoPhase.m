function S=UpstreamExpTwoPhase(Grid,S,Fluid1,Fluid0,V,q,T)

Nx=Grid.Nx; Ny=Grid.Ny; Nz=Grid.Nz;          % number of grid points
N=Nx*Ny*Nz;                                  % number of unknowns
pv = Grid.V(:).*Grid.por(:);                 % pore volume=cell volume*porosity

fi=max(q,0);                                 % inflow from wells
XP=max(V.x,0); XN=min(V.x,0);                % influx and outflux, x-faces
YP=max(V.y,0); YN=min(V.y,0);                % influx and outflux, y-faces
ZP=max(V.z,0); ZN=min(V.z,0);                % influx and outflux, z-faces

Vi = XP(1:Nx,:,:)+YP(:,1:Ny,:)+ZP(:,:,1:Nz)-...    % total flux into
XN(2:Nx+1,:,:)-YN(:,2:Ny+1,:)-ZN(:,:,2:Nz+1);      % each gridblock
pm = min(pv./(Vi(:)+fi));                    % estimate of influx
cfl = ((1-Fluid1.Sr-Fluid0.Sr)/3)*pm;        % CFL restriction
Nts = ceil(T/cfl);                           % number of local time steps
dtx = (T/Nts)./pv;                           % local time steps

A=GenSatMatrix(Grid,V,q);                    % system matrix - viscous sat. independent term
A=spdiags(dtx,0,N,N)*A;                      % A * dt/|Omega i|

fi=max(q,0).*dtx;                            % injection

fprintf('Number of local time-steps in UpstreamExpTwoPhase = %d\n', Nts);
for t=1:Nts
   [m1,m0]=RelativePerm(S,Fluid1,Fluid0);    % compute mobilities
   f1 = m1./(m1+m0);                         % compute fractional flow
   S = S+(A*f1+fi);
end