function C=Upstream(Grid,C,V,q,T)

Nx=Grid.Nx; Ny=Grid.Ny; Nz=Grid.Nz;          % number of grid points
N=Nx*Ny*Nz;                                  % number of unknowns

% mobile pore volume = cell volume*porosity*saturation 
pv = Grid.V(:).*Grid.por(:).*Grid.sat(:);    

fi=max(q,0);                                 % inflow from wells
XP=max(V.x,0); XN=min(V.x,0);                % influx and outflux, x-faces
YP=max(V.y,0); YN=min(V.y,0);                % influx and outflux, y-faces
ZP=max(V.z,0); ZN=min(V.z,0);                % influx and outflux, z-faces

Vi = XP(1:Nx,:,:)+YP(:,1:Ny,:)+ZP(:,:,1:Nz)-...    % total flux into
XN(2:Nx+1,:,:)-YN(:,2:Ny+1,:)-ZN(:,:,2:Nz+1);      % each gridblock
pm = min(pv./(Vi(:)+fi));                    % estimate of influx
cfl = (1/3)*pm;                              % CFL restriction
Nts = ceil(T/cfl);                           % number of local time steps
dtx = (T/Nts)./pv;                           % local time steps

A=GenA(Grid,V,q);                            % system matrix - convective term
A=spdiags(dtx,0,N,N)*A;                      % A * dt/|Omega i| 

fi=max(q,0).*dtx;                            % injection

fprintf('Number of local time-steps = %d\n', Nts);
for t=1:Nts
   % Concentration update for each local time step
   C = C + (A*C+fi);
end
