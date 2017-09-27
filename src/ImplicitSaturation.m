function [S,output] = ImplicitSaturation(Grid,S,Fluid1,Fluid0,V,q,T,opt)
% Fully implicit solve of the saturation equation for a given time step 
% using an adaptive time-stepping algorith combined with Newton-Raphson 
% method for inner nonlinear iterations.
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
%    opt.tol               - absolute tolerance for Newton-Raphson
%                            iterations as a convergence criterion
%    opt.maxit             - maximum number of allowed NR iterations 
%    opt.min_dt            - minimum allowed time step size below which the
%                            nonlinear NR iterations are stopped
%
% OUTPUTS:
% S                 - Array of cell-centered non-wetting phase saturation 
% output            - Output structure holding statistics on performance of
%                     the solver:
%    output.backsteps      - number of time backsteps 
%    output.nr_iterations  - total number of Newton-Raphson iterations 
%    output.dt             - selected time step by the adaptive algorithm 
%    output.convergence    - logical indicating the convergence status on
%                            exit (true or false)
%
% Author: M.A. Sbai, Ph.D.
%         BRGM (French Geological Survey) 
%         D3E  (Direction Eau, Environnement, Echotechnologies)
% 

% number of nonlinear system unknowns
N = Grid.Nx*Grid.Ny*Grid.Nz; 

% assemble advective system matrix
A = UpwindAdvectionMatrix(Grid,V,q);

% extract Newton-Raphson optional parametres 
if (nargin > 7)
   tol     = opt.tol;
   maxiter = opt.maxiter;
   min_dt  = opt.min_dt;
else 
   tol     = 1e-4;
   maxiter = 50;
   min_dt  = 1;
end

output.convergence = false; 
output.backsteps   = 0;

% savegard current saturations for later retreival if we do backstepping
% later
S00 = S;

%--- Main newton-raphson loop 
while ~output.convergence
    
    % current time step.
    dt  = T/2^(output.backsteps);
    
    % check that it is below the allowed minimum 
    if dt<min_dt 
        fprintf('Error: current time step is below the minimum allowed!\n');
        break;
    end
    
    % save current local time step in output structure 
    output.dt = dt;
    
    % time step divided by pore volume 
    dtx = dt./(Grid.V(:)*Grid.por(:));
    
    % extract 'dimensionless' injection flow rate
    fi = max(q,0).*dtx;
    
    % 'dimensionless' advection matrix
    B = spdiags(dtx,0,N,N)*A;

    % initialize counter of total Newton-Raphson iterations for all sub 
    % time steps
    output.nr_iterations = 0;
    
    % index of current sub-time step
    I = 0;
    
    %--- Sub-time steps loop
    while I < 2^(output.backsteps)                       
       
    % savegard saturations of previous time step 
    S0 = S; 
    
    % initialize the norm of saturations residuals 
    ds_norm = 1;
    
    % intilize Newton-Raphson iterations counter for this time sub-step
    nr_iterations = 0; 
    
    % increment local time step counter
    I = I+1;

        %--- Newton Raphson Iteration -------------------------------------
        while ds_norm > tol && nr_iterations < maxiter 
            
            % calculate fluids mobilities and their derivatives
            [M1,M0,dM1,dM0] = RelativePerm(S,Fluid1,Fluid0); 
            
            % calculate derivative of the fractional flow versus saturation
            % df/ds where f is the nonlinear fractional function
            df = dM1./(M1 + M0) - M1./(M1+M0).^2.*(dM1+dM0);
            
            % calculate dF(S) where F(S)=0 is the sparse nonlinear system
            % of algebraic equations being solved
            dF = speye(N) - B*spdiags(df,0,N,N);
            
            % calculate fractional flow term
            f = M1./(M1+M0);
            
            % calculate current residuals: F(S)
            F = S - S0 - (B*f + fi);
            
            % solve a linear system to get the NR increment dS
            dS = -dF\F;
            
            % update current iterate S
            S = S + dS;
            
            % calculate norm of increment dS 
            ds_norm = norm(dS);
            
            % update local NR iterations count
            nr_iterations = nr_iterations + 1;                                   
        
        end
        %------------------------------------------------------------------

        % check for convergence. If convergence if not yet attained
        % (residuals norm still below the prescribed tolerance) we do
        % backstepping (i.e. the time step solve is repeated with a smaller
        % time step size)
        if ds_norm > tol
            I = 2^(output.backsteps); 
            S = S00; 
        end
        
        % update total number of nonlinear iterations count
        output.nr_iterations = output.nr_iterations + nr_iterations;
        
    end
    
    % check for convergence.
    if ds_norm < tol 
        
        % successful exit of the adaptive algorithm
        output.convergence = true;
        
        if output.backsteps>0, fprintf('\n'); end
        fprintf('Converged in %d time backsteps and %d Newton-Raphson iterations\n', ...
                output.backsteps, output.nr_iterations);
   
    else
        
        % so, if not converged yet, decrease the local time step size
        output.backsteps = output.backsteps + 1;
        fprintf('.');
      
   end
   
end % adaptive local time step loop

end 
