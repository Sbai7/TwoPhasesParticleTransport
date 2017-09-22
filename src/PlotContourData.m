function PlotContourData(G, data, options, varargin)
% Plots iso-contours of cell-centered scalar data on a given 2D/3D cartesian 
% grid  
%
% INPUTS:
% G                 - Cartesian grid structure as created by calling 
%                     cartGrid routine   
% data              - array of scalar data to visualize on external facets
%                     of the grid 
% varargin (optional)
%                   - optional property-value pair arguments submitted to 
%                     contourf routine. Type 'help contourf' in the command 
%                     window to see the list of supported options
%
% Author: M.A. Sbai, Ph.D.
%         BRGM (French Geological Survey) 
%         D3E  (Direction Eau, Environnement, Echotechnologies)
%

data = data(:);
assert(G.ne==size(data,1), ...
    'Number of rows in data must equal the number of cells in the grid!');

ni = G.ni; nj = G.nj; 
Nx = ni-1; Ny = nj-1;
Dx = G.coord(ni,1)-G.coord(1,1);
Dy = G.coord(ni*nj,2)-G.coord(1,2);

if strcmp(options.fill,'on') 
    contourf(linspace((Dx/Nx)/2,Dx-(Dx/Nx)/2,Nx),...
             linspace((Dy/Ny)/2,Dy-(Dy/Ny)/2,Ny),...
             reshape(data,Nx,Ny)',varargin{:});
else
    contour( linspace((Dx/Nx)/2,Dx-(Dx/Nx)/2,Nx),...
             linspace((Dy/Ny)/2,Dy-(Dy/Ny)/2,Ny),...
             reshape(data,Nx,Ny)',varargin{:});
    
end
end

