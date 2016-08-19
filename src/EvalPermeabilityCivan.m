function k_factor=EvalPermeabilityCivan(phi0,phi,kappa,f,n)
%-------------------------------------------------------------------------------
% Calculates the reduced factors of permeability damage according to 
% Liu & Civan equation.
%
% phi0 is the initial gridded porosity
% phi is the current gridded porosity
% kappa is the residual permeability of plugged formation 
% f is the flow efficiency factor that is calculated from mass of particles that
%    are trapped at the pore throats. Should be calculated by the calling routine
% n is the exponent in Civan's law
%
% Author: M.A. Sbai, Ph.D.
% Creation date: 20/11/2008 
% Last revised:  20/11/2008 
%-------------------------------------------------------------------------------

% initialisation
gs = size(phi,1);         % size of gridded vectors
k_factor = zeros(gs,1);

% Calculation step
k_factor = (ones(gs,1)-f).*kappa + f.*(phi./phi0);

for i=1:gs
   k_factor(i,1) = k_factor(i,1)^n;
end
