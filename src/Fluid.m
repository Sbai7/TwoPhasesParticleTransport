classdef Fluid < handle
% A class for a fluid phase moving in a porous media. 

    properties (Access = public)
        name@char           % string name of this fluid (H2O, cO2, etc)
        density             % fluid phase density 
        viscosity           % fluid phase viscosity 
        Sr                  % residual fluid saturation
        compressibility     % fluid phase compressibility
    end
    
    
    methods (Access = public)
        
        function this = Fluid(name,in) 
        % Constructor function 
       
        assert(isa(name,'char'));
        this.name = name;
        
        in = in(:);
        assert(isa(in,'double')); 
        
        switch length(in)
   
        case 2       % initialise from density and viscosity only
            this.density   = in(1);
            this.viscosity = in(2);

        case 3       % initialise from density, viscosity, and Sr
            this.density   = in(1);
            this.viscosity = in(2);
            this.Sr        = in(3);

        case 4       % all member fields
            this.density   = in(1);
            this.viscosity = in(2);
            this.Sr        = in(3);
            this.compressibility = in(4);

        otherwise
            error('Wrong Fluid constructor.');
        
        end
        
        end
        
        %------------------------------------------------------------------
   
        % Other methods to implement in this class 
        % get_density(T,P) or get_density(T,P,S) 
        % get_viscosity(T,P) or get_viscosity(T,P,S)
        % get_compressibility(...)
        
        
    end
    
end
