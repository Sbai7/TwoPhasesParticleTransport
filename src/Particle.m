classdef Particle < handle
% A class for solid moving particles in porous media which are subject to 
% mechanistic processes of deposition (in pore bodies and throats), 
% hydrodynamic and colloidal release from pore bodies, and mass transfer 
% between fluid pahses of different mobilities. 
% 
% Author: M.A. Sbai
%         BRGM (French Geological Survey) 
%         D3E  (Direction Eau, Environnement, Echotechnologies)
% 

    properties (Access = public)
        diameter        % mean-size diamater 
        density         % density in kg/m^3 
        diff_coeff      % diffusion coefficient 
        uc              % critical fluid phase velocity 
        Csc             % critical fluid salinity 
        alpha_h         % hydrodynamic release rate 
        alpha_cl        % colloidal mobilisation rate 
        alpha_d         % deposition rate in pore surfaces (bodies) 
        alpha_pt        % deposition rate in pore throats 
        alpha_fe        % flow efficiency factor for civan's k-phi model 
        fluid@Fluid     % fluid phase to which the particles belong 
        N               % Dimension of arrays dependent on grid size
                        % like concentration, deposited conc's, etc ... 
        C = [];         % Initial mobile concentration
        C_dep = [];     % Initial deposited concentration on pore bodies
        C_pt  = [];     % Initial deposited concentration on pore throats
        C_tsf = [];     % Initial mass tranfer between two-phases
    end
    
    
    methods (Access = public)

        function this = Particle(fluid,param)
        % Constructor function 

        assert(isa(fluid,'Fluid'));
        this.fluid = fluid;
        
        param = param(:);
        assert(isa(param,'double'));
        assert(length(param)==10);
          
     	this.diameter   = param(1);
    	this.density    = param(2); 
      	this.diff_coeff = param(3);
      	this.uc         = param(4);
       	this.Csc        = param(5);
     	this.alpha_h    = param(6); 
      	this.alpha_cl   = param(7); 
      	this.alpha_d    = param(8); 
      	this.alpha_pt   = param(9); 
      	this.alpha_fe   = param(10);
   
        end

        %------------------------------------------------------------------

        function out = diffusion_coef(this,Tdegc,viscosity)
        % Returns the diffusion coefficient of the particle based on its 
        % geometric properties and temperature.
        % This models dilute suspensions of spheres as given by the 
        % well-known Stockes-Einstein equation.
        
        assert(isa(this,'Particle'));
        assert(isa(Tdegc,'double') && isscalar(Tdegc));
        assert(isa(viscosity,'double') && isscalar(viscosity));

        k = 1.3806504e-23;      % Boltzmann's constant (CODATA value) 
        TK = Tdegc - 273.15;	% conversion from degC to degK units  
        out = (k*TK)/(3*pi*viscosity*this.diameter);
        
        end

        %------------------------------------------------------------------

        function out = eff_diffusion_coef(this)
        %EFF_DIFFUSION_COEF returns the effective diffusion coefficient of  
        % the particle based on its geometric properties and temperature.
        % This models the effective diffusion coefficient for uniform or  
        % non-uniform sphere packing as published by Weissberg formula.  

        dim = size(phi,1);
        out = (diffusion(T,visc).*phi) ./ (ones(dim,1) - 0.5.*log(phi));

        end

        %------------------------------------------------------------------
        
        function [this,out] = transport(this,Grid,V,q,Cs,dt)
        % calculate the transport by convection-diffusion-kinetics 
        % of this particle over time step dt and updates concentrations 
        % of mobile particles and in pore surface/throat sites. 

        %---
        % solve the convective part of transport for this time-step.
        this.C = ImplicitConcentration(Grid,this.C,V,q,dt);
   
        %---
        % Solve for particle kinetics including the following mechanisms:
        % (1) Hydrodynamic release effect, (2) Colloidal release effect, 
        % (3) Surface deposition, (4) Pore-throat deposition, (5) interphase 
        % transfer.   
   
        % Update deposited mass of the particle & get their rates
        [this,diff1] = particleDeposits(Grid,V,this,Cs,dt);
   
        % Update pore-throat mass of the particle & get their rates
        [this,diff2]  = particlePtDeposits(Grid,V,this,dt);
   
        % Mass transfer between phases
        %   [C_tsf,diff3] = particleTransfer(Grid,V,this,C,C_tsf,dt);
        %   [this,diff3]  = particleTransfer(this.Grid,V,dt);
   
        % Update particle concentration from newly calculated total rate of 
        % mass transfer
        this.C = this.C - diff1 - diff2;
      
        %---
        % Solve the diffusive part of transport for this time-step 
        % --> reuse the TPFA solver (in transient mode) for this time step.
        [this,~] = this.diffusion(Grid,dt);
   
        out = 1;
        
        end
        
        %------------------------------------------------------------------
        
        function [this,Flux] = diffusion(this,Grid,dt)

            % Isotropic diffusion tensor 
            D = this.diff_coeff.*ones(3,Grid.Nx,Grid.Ny,Grid.Nz);

            % Compute particle concentration and extract diffusive fluxes
            [this.C,Flux] = Diffusion(Grid,D,this.C,dt);

        end
        
        %------------------------------------------------------------------
        
    end
    
end