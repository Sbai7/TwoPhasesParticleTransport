function S=ImplicitSaturation(Grid,S,Fluid1,Fluid0,V,q,T,opt);

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

conv=0; IT=0; S00=S;
backStep = 0;
while conv==0;
   dt  = T/2^IT; 			               % timestep
   dtx = dt./(Grid.V(:)*Grid.por(:)); 	% timestep / pore volume
   fi  = max(q,0).*dtx; 			      % injection
   B   = spdiags(dtx,0,N,N)*A;

   I=0;
   while I<2^IT;                       % loop over sub-timesteps
   S0=S; dsn=1; it=0; I=I+1;

      while dsn>tol & it<maxiter;  		% I T E R A T I O N
         [Mw,Mo,dMw,dMo]=RelativePerm(S,Fluid1,Fluid0);     % mobilities and derivatives
         df=dMw./(Mw + Mo)-Mw./(Mw+Mo).^2.*(dMw+dMo); % df w/ds
         dG=speye(N)-B*spdiags(df,0,N,N);             % G'(S)
         fw = Mw./(Mw+Mo);                            % fractional flow
         G = S-S0-(B*fw+fi);                          % G(s)
         ds = -dG\G;                                  % increment ds
         S = S+ds;                                    % update S
         dsn = norm(ds);                              % norm of increment
         it = it+1;                                   % number of N-R iterations
      end

      if dsn>tol; I=2^IT; S=S00; end % check for convergence
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
