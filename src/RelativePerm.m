function [M1,M0,dM1,dM0]=RelativePerm(s,Fluid1,Fluid0)

S = (s-Fluid1.Sr)/(1-Fluid1.Sr-Fluid0.Sr);   % Rescale saturations
M1 = S.^2/Fluid1.viscosity;                  % Fluid1 mobility
M0 =(1-S).^2/Fluid0.viscosity;               % Fluid0 mobility

if (nargout==4)
   dM1 = 2*S/Fluid1.viscosity/(1-Fluid1.Sr-Fluid0.Sr);
   dM0 = -2*(1-S)/Fluid0.viscosity/(1-Fluid1.Sr-Fluid0.Sr);
end