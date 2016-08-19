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
N = Grid.Nx*Grid.Ny*Grid.Nz;
Grid.N = N;

Grid.V =(Dx/Grid.Nx)*(Dy/Grid.Ny)*(Dz/Grid.Nz); % Cell volumes

% Gridded **scalar** permeability dataset
Grid.K = 0.85e-12.*ones(3,Grid.Nx,Grid.Ny,Grid.Nz);    % Unit-Darcy permeability
%%Grid.K = 10.^(-13*rand(3,Grid.Nx,Grid.Ny,Grid.Nz));

% Gridded porosity data
Grid.por = 0.3.*ones(Grid.Nx,Grid.Ny,Grid.Nz);

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
              0.03);               % Residual saturation

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
S = co2.Sr*ones(N,1);                       % Initial CO2 saturation
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

% setup NR options
opt.tol     = 1e-5;
opt.maxiter = 15;

in = index+10;
Pw=[co2.Sr; pt1.C(in); pt1.C_dep(in); pt1.C_pt(in); 1]; 
Pc=[co2.Sr; pt1.C(N); pt1.C_dep(N); pt1.C_pt(N)]; Tt=0;


%==========================================================================
% Begin time step loop (Now we assume same time stepping for flow &
% particle transport equations !)
%==========================================================================
tic
for t=1:nt
    
   fprintf('Solving two-phase flow problem. Time = %f days\n', t*dt/day);

   [P,V]=TwoPhasePressure(Grid,S,co2,water,Qw,P,dt);         % pressure solver
   
   %=======================================================================
   % Solve and plot CO2 saturation
   %=======================================================================
   %%S = UpstreamExpTwoPhase(Grid,S,co2,water,V,Qs,dt);
   S = ImplicitSaturation(Grid,S,co2,water,V,Qs,dt,opt);
   Grid.sat = reshape(S,Grid.Nx,Grid.Ny,Grid.Nz);
   
   %=======================================================================
   % Transport all particles in the system.
   %=======================================================================
   [m1,m0]=RelativePerm(S,co2,water);    
   [P1,V1]=TwoPhasePressure(Grid,m1.*S,co2,water,Qw,P,dt); 
   pt1.transport(Grid,V1,Qw,Cs,dt);
   
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
      
   % monitor injection well ... 
   Tt=[Tt, t*dt/day];
   Pw=[Pw, [S(in); pt1.C(in); pt1.C_dep(in); pt1.C_pt(in); k_ratio(in)]];
   subplot(1,2,1);
   plot(Tt, Pw(1,:), Tt, Pw(2,:), Tt, Pw(3,:), Tt, Pw(4,:), Tt, Pw(5,:));
   axis([0,dt*nt/day,-0.05,1.1]);
   legend('CO2 Saturation','Mobile conc','Surface conc','Pore throat conc','k/k0');
   drawnow;

   % monitor one production well ...
   Pc=[Pc, [S(N); pt1.C(N); pt1.C_dep(N); pt1.C_pt(N)]];
   subplot(1,2,2);
   plot(Tt, Pc(1,:), Tt, Pc(2,:), Tt, Pc(3,:), Tt, Pc(4,:));
   axis([0,dt*nt/day,-0.05,1.1]);
   legend('CO2 Saturation','Mobile conc','Surface conc','Pore throat conc');
   drawnow;
   
end

fprintf('Total simulation loop '); 
toc 
