classdef RecDomain < handle
% A class for simple description of a 3D rectangular domain.

    properties (Access = public)
        Lx;          % Rectangular domain extent along X-direction 
        Ly;          % Rectangular domain extent along Y-direction 
        Lz;          % Rectangular domain extent along Z-direction 
    end

   
    methods (Access = public)
       
        function this = RecDomain(L)
        % class constructor: simply creates a domain with given parameters.

        L = L(:);
        assert(isa(L,'double'), ...
            'Input argument must be a scalar or 3 elements double array.');
        
        switch length(L) 
            case 1
                this.Lx = L; 
                this.Ly = L; 
                this.Lz = L; 
      
            case 3
                this.Lx = L(1);
                this.Ly = L(2);
                this.Lz = L(3);

            otherwise 
                error('Wrong RecDomain constructor!');
        end
        
        end
        
        %------------------------------------------------------------------

   end
   
end