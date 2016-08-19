function k_factor=EvalPermeabilityWennberg(phi0,phi,beta,xsi,d_pores,d_grains)
%-------------------------------------------------------------------------------
% Calculates the reduced factors of permeability damage according to 
% Wennberg & Sharma equation. This equation splits the permeability damage to 
% three factors.
%
% phi0 is the initial gridded porosity
% phi is the current gridded porosity
% beta is the 'distributed' damage factor of porous materials. We allow a vector
%    entry form for this parameter such that materials with different damage 
%    factors could be modelled easily by this routine. 
% xsi is the volume fraction of deposited particles
% d_pores is the mean diameter of pores 
% d_grains is the mean diameter of 'spherical' grains of porous materials.
%
% Author: M.A. Sbai, Ph.D.
% Creation date: 20/11/2008 
% Last revised:  20/11/2008 
%-------------------------------------------------------------------------------


% initialisation ...
gs = size(phi,1);         % size of gridded vectors
k_factor = zeros(gs,1);

% Calculate contribution from porosity change 
k_factor = (phi.^3).*((ones(gs,1)-phi0).^2) ...
           ./(phi0.^3)./((ones(gs,1)-phi).^2);
           
% Calculate contribution from tortuosity change
k_factor = k_factor./(ones(gs,1)+beta.*xsi);

% Calculate contribution from the change in specific surface area. 
k_factor = k_factor.*(ones(gs,1)+(xsi./(ones(gs,1)-phi0))) ...
           ./(ones(gs,1)+(xsi.*d_grains./(ones(gs,1)-phi0)./d_pores));
           