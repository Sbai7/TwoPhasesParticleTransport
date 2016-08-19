function [C,Flux]=diffusion(Grid,particle,C,dt)

% Isotropic diffusion tensor 
D = particle.diff_coeff.*ones(3,Grid.Nx,Grid.Ny,Grid.Nz);

% Compute particle concentration and extract diffusive fluxes
[C,Flux] = fv_diffusion(Grid,D,C,dt);
