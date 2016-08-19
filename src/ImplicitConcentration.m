function C=ImplicitConcentration(Grid,C,V,q,T,opt);

N = Grid.Nx*Grid.Ny*Grid.Nz; 		% number of unknowns
A = GenSatMatrix(Grid,V,q); 	   % system matrix

% extract NR options 
if (nargin > 7)
   tol     = opt.tol;
   maxiter = opt.maxiter;
else 
   tol     = 1e-3;
   maxiter = 10;
end

conv=0; IT=0; C00=C;
backStep = 0;
while conv==0;
   dt  = T/2^IT; 			               % timestep
   dtx = dt./(Grid.V(:).*Grid.por(:).*Grid.sat(:)); 	% timestep / pore volume
   fi  = max(q,0).*dtx; 			      % injection
   B   = spdiags(dtx,0,N,N)*A;

   I=0;
   while I<2^IT;                       % loop over sub-timesteps
   C0=C; dsn=1; it=0; I=I+1;

      while dsn>tol & it<maxiter;  		% I T E R A T I O N
         df=ones(N,1);                                % df/dc = 1
         dG=speye(N)-B*spdiags(df,0,N,N);             % G'(S)
         G = C-C0-(B*C+fi);                           % G(s)
         dc = -dG\G;                                  % increment ds
         C = C+dc;                                    % update S
         dsn = norm(dc);                              % norm of increment
         it = it+1;                                   % number of N-R iterations
      end

      if dsn>tol; I=2^IT; C=C00; end % check for convergence
   end

   if dsn < tol;                     % check for convergence
      conv = 1;
      if IT>0; fprintf('\n'); end
      fprintf('Newton-Raphson iteration converged in %d steps\n', IT+1);
   else                              % if not converged, decrease sub-timestep
      IT = IT + 1;
      if backStep==1;
         fprintf('.');
      else
         fprintf('back-stepping .');
      end
      backStep = 1;
   end
end % timestep decrease by factor 2
