clear all; % Useful to clear stack of all previous variables

% Domain size along X, Y & Z directions
dom = RecDomain(365.76,670.56,0.6096);
Dx = dom.Dx;
Dy = dom.Dy;
Dz = dom.Dz;

% Number of cells along X, Y & Z directions
Grid.Nx = 60;
Grid.Ny = 220;
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

% Read Gridded permeability map
Grid.K = 1.4e-11.*ones(3,Grid.Nx,Grid.Ny,Grid.Nz);    % Unit-Darcy permeability
%%load 'spe_layer1_perm.mat' perm_layer;
%%Grid.K = perm_layer*9.869232667160130e-13;
%%kk = Grid.K(:,1);
%%Grid.K = reshape(Grid.K',3,Nx,Ny,Nz); 
   
% Read Gridded porosity map
%%Grid.por = 0.3.*ones(Grid.Nx,Grid.Ny,Grid.Nz);
load 'spe_layer1_phi.mat' phi_layer;
phi = phi_layer;
phi = max (phi(:), 1e-3); % to avoid zero porosity cells
Grid.por = reshape(phi',Nx,Ny,Nz); 
   
figure 
   subplot(3,2,1);
   % plot filled contours at the midpoints of the grid cells
   contourf(linspace((Dx/Grid.Nx)/2,Dx-(Dx/Grid.Nx)/2,Grid.Nx),...
            linspace((Dy/Grid.Ny)/2,Dy-(Dy/Grid.Ny)/2,Grid.Ny),...
            reshape(phi,Nx,Ny)',21);

   axis('equal'); caxis('auto');       % equal axes and color
   colorbar; title('(-)');
   drawnow;                            % force update of plot

%%clear spe_layer36_perm kk;
clear spe_layer36_phi;

% Gridded 'total' compressibility data
Grid.compr = 1e-5.*ones(Grid.Nx,Grid.Ny,Grid.Nz);

% Cell-centered injection/production flow rates
Total      = 150/3600;  % m^3/h >> m^3/s
Qw         = zeros(N,1);
Qw([1 N])  = [-Total/4 -Total/4];
Qw(Nx)     = -Total/4;
Qw(N-Nx+1) = -Total/4;
index      = (Nx*(Ny-1))/2;
Qw(index)  = Total;

Total      = 150/3600;  % m^3/h >> m^3/s
Qs         = zeros(N,1);
Qs([1 N])  = [-Total/4 -Total/4];
Qs(Nx)     = -Total/4;
Qs(N-Nx+1) = -Total/4;
index      = (Nx*(Ny-1))/2;
Qs(index)  = Total;

% Fluid properties
oil        = Fluid(1200.0, ...          % Density en kg/m^3
                   3.e-3, ...           % Viscosity en Pa.s
                   0.3);                % Residual saturation
co2        = Fluid(870.2362, ...        % Density en kg/m^3
                   0.0768e-3, ...       % Viscosity en Pa.s
                   0.05);               % Residual saturation

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

% Initialisation of concentration arrays => particle 1
pt1.C  = zeros(N,1);               % Initial concentration
pt1.C_dep = zeros(N,1);            % Initial surf. deposited conc
pt1.C_pt  = zeros(N,1);            % Initial throats deposited conc
pt1.C_tsf = zeros(N,1);            % Initial mass tranfer <=> phases
Cs = zeros(N,1);                   % Initial salt concentration

% Initialisation of the pressure array
P = 300e5.*ones(Grid.Nx,Grid.Ny,Grid.Nz);   % Initial fluid pressure
S = 0.3*ones(N,1);                          % Initial CO2 saturation
Grid.sat = reshape(S,Grid.Nx,Grid.Ny,Grid.Nz); % Gridded fluid saturation


% Time stepping, etc ...
day = 3600*24;                     % seconds/day
nt  = 30;                          % Time steps
dt  = 30*day/nt;                   % ...

% Copy of initial porosity & permeability
Grid.por0 = Grid.por;
Grid.K0   = Grid.K;
% Other parameters
beta      = 0.*ones(N,1);          % distributed damage factor (Wennberg k-phi model)
kappa     = 0;                     % Residual K accounting for the valve effect (Civan k-phi model)
d_pores   = 3e-5; 
d_grains  = 0.175e-3; 

%==========================================================================
% Begin time step loop (Now we assume same time stepping for flow &
% particle transport equations !)
%==========================================================================
tic
for t=1:nt
    
   fprintf('Solving two-phase flow problem. Time = %f days\n', t*dt/day);

      [P,V]=TwoPhasePressure(Grid,S,co2,oil,Qw,P,dt);         % pressure solver
%%   subplot(3,2,1);
   
   % plot pressure filled contours at the midpoints of the grid cells
%%   contourf(linspace((Dx/Grid.Nx)/2,Dx-(Dx/Grid.Nx)/2,Grid.Nx),...
%%            linspace((Dy/Grid.Ny)/2,Dy-(Dy/Grid.Ny)/2,Grid.Ny),...
%%            reshape(P/1e5,Grid.Nx,Grid.Ny),11);

%%   axis('equal'); caxis('auto');       % equal axes and color
%%   colorbar; title('(A)');
%%   drawnow;                            % force update of plot
   
   %=======================================================================
   % Solve and plot CO2 saturation
   %=======================================================================
   S = UpstreamExpTwoPhase(Grid,S,co2,oil,V,Qs,dt);
   Grid.sat = reshape(S,Grid.Nx,Grid.Ny,Grid.Nz);
   
   %=======================================================================
   % Transport all particles in the system.
   %=======================================================================
   [m1,m0]=RelativePerm(S,co2,oil); 
   %%f1 = m1./(m1+m0); f1 = reshape(f1,Grid.Nx,Grid.Ny,Grid.Nz);
   %%V1.x = zeros(Nx+1,Ny,Nz);
   %%V1.y = zeros(Nx,Ny+1,Nz);
   %%V1.z = zeros(Nx,Ny,Nz+1);
   
   % Calculate Velocity of the CO2 fluid phase
   %%V1.x(2:Nx+1,:,:) = f1(:,:,:).*V.x(2:Nx+1,:,:); V1.x(1,:,:) = V.x(1,:,:);
   %%V1.y(:,2:Ny+1,:) = f1(:,:,:).*V.y(:,2:Ny+1,:); V1.y(:,1,:) = V.y(:,1,:);
   %%V1.z(:,:,2:Nz+1) = f1(:,:,:).*V.z(:,:,2:Nz+1); V1.z(:,:,1) = V.z(:,:,1);
   
   [P1,V1]=TwoPhasePressure(Grid,m1.*S,co2,oil,Qw,P,dt); 
   pt1.transport(Grid,co2,V1,Qw,Cs,dt);
   
   %=======================================================================
   % Evaluate porosity and permeability chnage.
   %=======================================================================
   Grid.por0 = reshape(Grid.por0,N,1);
   Grid.por  = reshape(Grid.por,N,1);
   Grid.por = Grid.por0 - ( (pt1.C_dep + pt1.C_pt) / pt1.density);
   
%   k_ratio = EvalPermeabilityWennberg(Grid.por0,Grid.por,beta, ...
%                                      (pt1.C_dep + pt1.C_pt)./pt1.C, ...
%                                      d_pores,d_grains);

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
   
   subplot(3,2,2);  
   % plot saturation filled contours at the midpoints of the grid cells
   contourf(linspace((Dx/Grid.Nx)/2,Dx-(Dx/Grid.Nx)/2,Grid.Nx),...
            linspace((Dy/Grid.Ny)/2,Dy-(Dy/Grid.Ny)/2,Grid.Ny),...
            reshape(S,Grid.Nx,Grid.Ny)',11);

   axis('equal'); caxis('auto');       % equal axes and color
   colorbar; title('(A)');
   drawnow;                            % force update of plot
   
   
   subplot(3,2,3);
   % plot filled contours at the midpoints of the grid cells
   contourf(linspace((Dx/Grid.Nx)/2,Dx-(Dx/Grid.Nx)/2,Grid.Nx),...
            linspace((Dy/Grid.Ny)/2,Dy-(Dy/Grid.Ny)/2,Grid.Ny),...
            reshape(pt1.C,Grid.Nx,Grid.Ny)',11);

   axis('equal'); caxis('auto');       % equal axes and color
   colorbar; title('(B)');
   drawnow;                            % force update of plot
   
   
   subplot(3,2,4);
   % plot filled contours at the midpoints of the grid cells
   contourf(linspace((Dx/Grid.Nx)/2,Dx-(Dx/Grid.Nx)/2,Grid.Nx),...
            linspace((Dy/Grid.Ny)/2,Dy-(Dy/Grid.Ny)/2,Grid.Ny),...
            reshape(pt1.C_dep,Grid.Nx,Grid.Ny)',11);

   axis('equal'); caxis('auto');       % equal axes and color
   colorbar; title('(C)');
   drawnow;                            % force update of plot
   
   
   subplot(3,2,5);
   % plot filled contours at the midpoints of the grid cells
   contourf(linspace((Dx/Grid.Nx)/2,Dx-(Dx/Grid.Nx)/2,Grid.Nx),...
            linspace((Dy/Grid.Ny)/2,Dy-(Dy/Grid.Ny)/2,Grid.Ny),...
            reshape(pt1.C_pt,Grid.Nx,Grid.Ny)',11);

   axis('equal'); caxis('auto');       % equal axes and color
   colorbar; title('(D)');
   drawnow;                            % force update of plot
   
   
   subplot(3,2,6);
   % plot filled contours at the midpoints of the grid cells
   contourf(linspace((Dx/Grid.Nx)/2,Dx-(Dx/Grid.Nx)/2,Grid.Nx),...
            linspace((Dy/Grid.Ny)/2,Dy-(Dy/Grid.Ny)/2,Grid.Ny),...
            reshape(k_ratio,Grid.Nx,Grid.Ny)',11);

   axis('equal'); caxis('auto');       % equal axes and color
   colorbar; title('(E)');
   drawnow;                            % force update of plot
   
   
 %%  subplot(3,2,6);
   % plot filled contours at the midpoints of the grid cells
 %%  contourf(linspace((Dx/Grid.Nx)/2,Dx-(Dx/Grid.Nx)/2,Grid.Nx),...
 %%           linspace((Dy/Grid.Ny)/2,Dy-(Dy/Grid.Ny)/2,Grid.Ny),...
 %%           reshape(Grid.por./Grid.por0,Grid.Nx,Grid.Ny),11);

 %%  axis('equal'); caxis('auto');       % equal axes and color
 %%  colorbar
 %%  drawnow;                            % force update of plot
   
end

fprintf('Total simulation loop '); 
toc 
