function S = ImplicitSaturation(Grid,S,Fluid1,Fluid0,V,q,T,opt)
% Fully implicit solve of the saturation equation for a given time step 
% using Newton-Raphson method.
%
% INPUTS:
% Grid              - Grid used for discretization 
% S                 - Array of cell-centered non-wetting phase saturation 
%                     at previous time step
% Fluid1            - Fluid properties of the non-wetting (injected) phase 
% Fluid0            - Fluid properties of the wetting (residant) phase
% V                 - Structure of mid-cells edges velocities along x,y,and
%                     z directions
% q                 - Injection/Production flow rates 
% T                 - Global time step 
% opt (optional)    - Newton-Raphson solver options
%
% OUTPUTS:
% S                 - Array of cell-centered non-wetting phase saturation 
%
% Author: M.A. Sbai, Ph.D.
%         BRGM (French Geological Survey) 
%         D3E  (Direction Eau, Environnement, Echotechnologies)
% 

% number of nonlinear system unknowns
N = Grid.Nx*Grid.Ny*Grid.Nz; 

% assemble convective system matrix
A = UpwindAdvectionMatrix(Grid,V,q);

% extract Newton-Raphson optional parametres 
if (nargin > 7)
   tol     = opt.tol;
   maxiter = opt.maxiter;
else 
   tol     = 1e-3;
   maxiter = 10;
end

converged = false; 
IT        = 0; 
S00       = S;

%--- Main newton-raphson loop 
while ~converged
    
    dt  = T/2^IT; 			               % current time step
    dtx = dt./(Grid.V(:)*Grid.por(:)); 	   % timestep / pore volume
    fi  = max(q,0).*dtx; 			       % injection
    B   = spdiags(dtx,0,N,N)*A;

    % total number of Newton-Raphson iterations for all time sub-steps
    total_nr_iterations = 0;
    
    I = 0;
    while I < 2^IT                       % loop over sub-timesteps
    
    S0              = S; 
    dx_norm         = 1;
    
    % number of Newton-Raphson iterations for this time sub-step
    nr_iterations   = 0; 
    I               = I+1;

        %--- Newton Raphson Iteration -------------------------------------
        while dx_norm > tol && nr_iterations < maxiter 
            
            % calculate fluids mobilities and their derivatives
            [M1,M0,dM1,dM0] = RelativePerm(S,Fluid1,Fluid0); 
            
            % calculate df/ds where f is the fractional flow term
            df = dM1./(M1 + M0) - M1./(M1+M0).^2.*(dM1+dM0);
            
            % calculate dF(S) where F is the zero sparse function
            dF = speye(N) - B*spdiags(df,0,N,N);
            
            % calculate fractional flow term
            f = M1./(M1+M0);
            
            % calculate current residuals: F(S)
            F = S - S0 - (B*f + fi);
            
            % solve sublinear system to get the saturation increment dS
            dS = -dF\F;
            
            % update current iterate S
            S = S + dS;
            
            % calculate norm of increment dS 
            dx_norm = norm(dS);
            
            % update nonlinear iterations count
            nr_iterations = nr_iterations + 1;                                   
        
        end
        %------------------------------------------------------------------

        % check for convergence
        if dx_norm > tol
            I = 2^IT; 
            S = S00; 
        end
        
        % update total number of nonlinear iterations count
        total_nr_iterations = total_nr_iterations + nr_iterations;
        
    end
    
    % check for convergence
    if dx_norm < tol   
        
        converged = true;
        if IT>0, fprintf('\n'); end
        fprintf('Converged in %d time sub-steps and %d Newton-Raphson iterations\n', ...
                IT+1, total_nr_iterations);
   
    else                          % if not converged, decrease sub-timestep
        
      IT = IT + 1;
      fprintf('.');
      
   end
   
end % timestep decrease by factor 2

end 
