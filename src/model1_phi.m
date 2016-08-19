clear all; % Useful to clear stack of all previous variables

% Domain size along X, Y & Z directions
dom = RecDomain(500,500,5);
Dx = dom.Dx;
Dy = dom.Dy;
Dz = dom.Dz;

% Number of cells along X, Y & Z directions
Grid.Nx = 101;
Grid.Ny = 101;
Grid.Nz = 1; 
Nx = Grid.Nx; Ny = Grid.Ny; Nz = Grid.Nz; 

% Calculate Grid cells sizes 
Grid.hx = (dom.Dx/Grid.Nx);
Grid.hy = (dom.Dy/Grid.Ny);
Grid.hz = (dom.Dz/Grid.Nz);

% Total number of grid blocks
N = Grid.Nx*Grid.Ny;
Grid.N = N;

Grid.V =(Dx/Grid.Nx)*(Dy/Grid.Ny)*(Dz/Grid.Nz); % Cell volumes

% Gridded 'total' compressibility data
Grid.compr = 1e-5.*ones(Grid.Nx,Grid.Ny,Grid.Nz);

% Cell-centered injection/production flow rates
Total = 150/3600;  % m^3/h >> m^3/s
Qw = zeros(N,1);
Qw([1 N])=[-Total/4 -Total/4];
Qw(Grid.Nx) = -Total/4;
Qw(N-Grid.Nx+1) = -Total/4;
index = ((Grid.Nx*Grid.Nx)-1)/2;
Qw(index) = Total;

Total = 150/3600;  % m^3/h >> m^3/s
Qs = zeros(N,1);
Qs([1 N])=[-Total/4 -Total/4];
Qs(Grid.Nx) = -Total/4;
Qs(N-Grid.Nx+1) = -Total/4;
index = ((Grid.Nx*Grid.Nx)-1)/2;
Qs(index) = Total;

% Fluid properties
water = Fluid(1065.5, ...          % Density en kg/m^3
              0.5065e-3, ...       % Viscosity en Pa.s
              0.3);                % Residual saturation
co2   = Fluid(870.2362, ...        % Density en kg/m^3
              0.0768e-3, ...       % Viscosity en Pa.s
              0.3);                % Residual saturation

% Particle properties 
pt1 = Particle(2e-6, ...           % Mean-size diameter => for particle population
               2500, ...           % Particles density
               1e-9, ...           % Diffusion coeff => important for small size particles
               0.2e-4, ...         % Critical fluid velocity
               0.1, ...            % Critical fluid salinity
               3.8e-4, ...         % Hydrodynamic release rate
               0, ...              % Colloidal mobilisation rate
               1.2e-4, ...         % Deposition rate in pore surfaces
               6.2e-6, ...         % Deposition rate in pore throats
               0.6, ...            % For civan k_phi model
               60, ...             % system temperature
               0.0768e-3 ...       % viscosity of the carying fluid phase
               );

% Time stepping, etc ...
day = 3600*24;                     % seconds/day
nt  = 30;                          % Time steps
dt  = 30*day/nt;                   % ...

% Other parameters
beta      = 0.*ones(N,1);          % distributed damage factor (Wennberg k-phi model)
kappa     = 0;                     % Residual K accounting for the valve effect (Civan k-phi model)
d_pores   = 3e-5; 
d_grains  = 0.175e-3; 

%==========================================================================
% Begin time step loop (Now we assume same time stepping for flow &
% particle transport equations !)
%==========================================================================
figure

tic
for ik=1:6

% Gridded permeability data
Grid.K = 0.85e-12.*ones(3,Grid.Nx,Grid.Ny,Grid.Nz);

% Initialisation of concentration arrays => particle 1
pt1.C  = zeros(N,1);               % Initial concentration
pt1.C_dep = zeros(N,1);            % Initial surf. deposited conc
pt1.C_pt  = zeros(N,1);            % Initial throats deposited conc
pt1.C_tsf = zeros(N,1);            % Initial mass tranfer <=> phases
Cs = zeros(N,1);                   % Initial salt concentration

% Initialisation of the pressure array
P = 300e5.*ones(Grid.Nx,Grid.Ny,Grid.Nz);   % Initial fluid pressure
S = 0.3*ones(N,1);                         % Initial CO2 saturation
Grid.sat = reshape(S,Grid.Nx,Grid.Ny,Grid.Nz); % Gridded fluid saturation

if ik==1
   Grid.por = 0.35.*ones(Grid.Nx,Grid.Ny,Grid.Nz);
elseif ik==2
   Grid.por = 0.30.*ones(Grid.Nx,Grid.Ny,Grid.Nz);
elseif ik==3
   Grid.por = 0.25.*ones(Grid.Nx,Grid.Ny,Grid.Nz);
elseif ik==4
   Grid.por = 0.20.*ones(Grid.Nx,Grid.Ny,Grid.Nz);
elseif ik==5
   Grid.por = 0.15.*ones(Grid.Nx,Grid.Ny,Grid.Nz);
else
   Grid.por = 0.10.*ones(Grid.Nx,Grid.Ny,Grid.Nz);
end                                 
  
% Copy of initial porosity & permeability
Grid.por0 = Grid.por;
Grid.K0   = Grid.K;
  
for t=1:nt
    
   fprintf('Solving two-phase flow problem. Time = %f days\n', t*dt/day);
   [P,V]=TwoPhasePressure(Grid,S,co2,water,Qw,P,dt);      % pressure solver
   
   %=======================================================================
   % Solve and plot CO2 saturation
   %=======================================================================
   S = UpstreamExpTwoPhase(Grid,S,co2,water,V,Qs,dt);
   Grid.sat = reshape(S,Grid.Nx,Grid.Ny,Grid.Nz);
   
   %=======================================================================
   % Transport all particles in the system.
   %=======================================================================
   [m1,m0]=RelativePerm(S,co2,water);
   [P1,V1]=TwoPhasePressure(Grid,m1.*S,co2,water,Qw,P,dt); 
   pt1.transport(Grid,co2,V1,Qw,Cs,dt);
   
   %=======================================================================
   % Evaluate porosity and permeability chnage.
   %=======================================================================
   Grid.por0 = reshape(Grid.por0,N,1);
   Grid.por  = reshape(Grid.por,N,1);
   Grid.por = Grid.por0 - ( (pt1.C_dep + pt1.C_pt) / pt1.density);

   f = abs(ones(N,1) - pt1.alpha_fe.*pt1.C_pt); % flow efficiency factor *** distributed ***                               
   k_ratio = EvalPermeabilityCivan(Grid.por0,Grid.por,kappa,f,3);
   
   Grid.por0 = reshape(Grid.por0,Grid.Nx,Grid.Ny,Grid.Nz);
   Grid.por  = reshape(Grid.por,Grid.Nx,Grid.Ny,Grid.Nz);
   
   k_ratio = reshape(k_ratio,Grid.Nx,Grid.Ny,Grid.Nz);
   for l=1:3
      for i=1:Grid.Nx
      for j=1:Grid.Ny
      for k=1:Grid.Nz
         Grid.K(l,i,j,k) = Grid.K0(l,i,j,k) * k_ratio(i,j,k);
      end
      end
      end
   end
   k_ratio = reshape(k_ratio,N,1);
  
   %=======================================================================
   % Plotting commands 
   %=======================================================================   
   subplot(3,2,ik);
   % plot filled contours at the midpoints of the grid cells
   contourf(linspace((Dx/Grid.Nx)/2,Dx-(Dx/Grid.Nx)/2,Grid.Nx),...
            linspace((Dy/Grid.Ny)/2,Dy-(Dy/Grid.Ny)/2,Grid.Ny),...
            reshape(k_ratio,Grid.Nx,Grid.Ny),11);

   axis('equal'); caxis('auto');       % equal axes and color
   colorbar; title('(E)');
   drawnow;                            % force update of plot
   
end
fprintf('End of section : %d \n\n',ik);
end

fprintf('Total simulation loop '); 
toc 
