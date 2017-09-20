function C = ImplicitConcentration(Grid,C,V,q,T,opt)
% Fully implicit solve of the 3D convective part of mass transport equation 
% for a given time step using Newton-Raphson method.
%
% INPUTS:
% Grid              - Grid used for discretization 
% C                 - Array of cell-centered concentration at previous 
%                     time step
% V                 - Structure of mid-cells edges velocities along x,y,and
%                     z directions
% q                 - Injection/Production flow rates 
% T                 - Global time step
% opt (optional)    - Newton-Raphson solver options
%
% OUTPUTS:
% S                 - Array of cell-centered concentration 
%
% Author: M.A. Sbai, Ph.D.
%         BRGM (French Geological Survey) 
%         D3E  (Direction Eau, Environnement, Echotechnologies)
% 

N = Grid.Nx*Grid.Ny*Grid.Nz;         % number of unknowns
A = UpwindAdvectionMatrix(Grid,V,q); % system matrix

% extract NR options 
if (nargin > 5)
   tol     = opt.tol;
   maxiter = opt.maxiter;
else 
   tol     = 1e-3;
   maxiter = 10;
end

conv=0; IT=0; C00=C;
while conv==0;
   dt  = T/2^IT; 			               % timestep
   dtx = dt./(Grid.V(:).*Grid.por(:).*Grid.sat(:)); 	% timestep / pore volume
   fi  = max(q,0).*dtx; 			      % injection
   B   = spdiags(dtx,0,N,N)*A;

   I=0;
   while I<2^IT;                       % loop over sub-timesteps
   C0=C; dsn=1; it=0; I=I+1;

      while dsn > tol && it < maxiter
         df=ones(N,1);                                % df/dc = 1
         dG=speye(N)-B*spdiags(df,0,N,N);             
         G = C-C0-(B*C+fi);                           
         dc = -dG\G;                                  
         C = C+dc;                                    
         dsn = norm(dc);                              
         it = it+1;                                   
      end

      if dsn>tol; I=2^IT; C=C00; end % check for convergence
   end

   if dsn < tol                   % check for convergence
      conv = 1;
      if IT>0; fprintf('\n'); end
      fprintf('Newton-Raphson iteration converged in %d steps\n', IT+1);
   else                           % if not converged, decrease sub-timestep
      IT = IT + 1;
      fprintf('.');
   end
end % timestep decrease by factor 2

end