function k_factor = EvalPermeabilityCivan(phi0,phi,kappa,f,n)
% Calculates the reduced factors of permeability damage according to 
% Liu & Civan (1993) correlation.
%
% INPUTS:
% phi0              - initial gridded porosity
% phi               - current gridded porosity
% kappa             - residual permeability of plugged formation 
% f                 - flow efficiency factor that is calculated from mass 
%                     of trapped particles at pore throats. This should be 
%                     calculated by the calling routine
% n                 - exponent in Civan's law
%
% Author: M.A. Sbai, Ph.D. 
%         BRGM (French Geological Survey) 
%         D3E  (Direction Eau, Environnement, Echotechnologies)
% 

% initialisation
gs = size(phi,1);         % size of gridded vectors

% Calculation step
k_factor = (ones(gs,1)-f).*kappa + f.*(phi./phi0);

for i=1:gs
   k_factor(i,1) = k_factor(i,1)^n;
end

end