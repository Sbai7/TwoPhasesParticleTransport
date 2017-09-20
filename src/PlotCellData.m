function PlotCellData(g, data, varargin)
% Plots cell-centered scalar data on a given 2D/3D cartesian grid  
%
% INPUTS:
% g                 - Cartesian grid structure as created by calling 
%                     cartGrid routine   
% data              - array of scalar data to visualize on external facets
%                     of the grid 
% varargin (optional)
%                   - optional property-value pair arguments submitted to 
%                     PlotGrid routine and subsequently to MATLAB patch 
%                     command. Type 'help patch' in the command window to 
%                     see the list of supported options
%
% Author: M.A. Sbai, Ph.D.
%         BRGM (French Geological Survey) 
%         D3E  (Direction Eau, Environnement, Echotechnologies)
% 

assert(g.ne==size(data,1), ...
    'Number of rows in data must equal the number of cells in the grid!');

% save data min/max
data_min = min(data); 
data_max = max(data); 

% transform scalar data to RGB values 
cmap = colormap();
rng  = [min(data), max(data)];
nc   = size(cmap, 1);

if ~(rng(2) > rng(1)),
    data = cmap(ceil(nc / 2), :);
else
	ix   = ceil(nc * (data(:) - rng(1)) ./ diff(rng));
	data = cmap(max(1, min(ix, nc)), :);
end

PlotGrid(g,'FaceColor','flat','CData',reshape(data,g.ne,1,3),varargin{:});

set(gca,'CLim',[data_min data_max]), colorbar

end