function [M1,M0,dM1,dM0] = RelativePerm(S_nw,Fluid1,Fluid0)
% Calculates fluids mobilities and their first-order derivatives. 
%
% INPUTS:
% S_nw              - Array of cell-centered fluid1 (non-wetting) phase 
%                     saturation 
% Fluid1            - Fluid properties of the non-wetting (injected) phase 
% Fluid0            - Fluid properties of the wetting (residant) phase
%
% OUTPUTS:
% M1                 - Array of cell-centered mobilities of Fluid1
% M0                 - Array of cell-centered mobilities of Fluid0
% dM1 (optional)     - Array of cell-centered derivatives of Fluid1
%                      mobility versus S_nw saturation (dM1/dS) 
% dM0 (optional)     - Array of cell-centered derivatives of Fluid0
%                      mobility versus S_nw saturation (dM0/dS)
%
% Author: M.A. Sbai, Ph.D.
%         BRGM (French Geological Survey) 
%         D3E  (Direction Eau, Environnement, Echotechnologies)
% 

%--- here we assume a relative permeability function kr(s) = s^2

S_nw = S_nw(:);

% scaled saturations
S_star = (S_nw-Fluid1.Sr)/(1-Fluid1.Sr-Fluid0.Sr);   

% mobility of Fluid1 
M1 = S_star.^2 ./ Fluid1.viscosity(:);          

% mobility of Fluid0
M0 =(ones(size(S_star))-S_star).^2 ./ Fluid0.viscosity(:); 

if (nargout==4)
   dM1 = 2*S_star./Fluid1.viscosity/(1-Fluid1.Sr-Fluid0.Sr);
   dM0 = -2*(ones(size(S_star))-S_star)./Fluid0.viscosity/ ... 
                                        (1-Fluid1.Sr-Fluid0.Sr);
end

end