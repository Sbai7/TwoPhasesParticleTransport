function G = cartGrid(x,y,z)
% Create a uniform/non-uniform cartesian grid from rectilinear arrays along
% each dimension.
%
% G = cartGrid(x)
% G = cartGrid(x,y)
% G = cartGrid(x,y,z) 
%
% Examples:
% 1) Create a 2D cartesian uniform grid over the domain 0 < x < 1000 and 0
% < y < 2000 with a 10 and 20m spacing, rspectively:
% >> G = cartGrid(0:10:1000, 0:20:2000); 
% >> plotGrid(G, 'FaceColor', 'none', 'EdgeColor', 'r');
% 
% 2) Create a 3D unifrom grid in X and Y, but non-uniform in Z over the
% domain 0 < x < 1000, 0 < y < 2000, and -50 < z < 25:
% >> G = cartGrid(0:10:1000, 0:20:2000, [25 20 15 10 5 0 -5 -15 -30 -50]);
% >> plotGrid(G, 'FaceColor', 'c', 'FaceAlpha', 0.7); 
% 
% Author: M.A. Sbai, Ph.D.
%         BRGM (French Geological Survey) 
%         D3E  (Direction Eau, Environnement, Echotechnologies)
% 

% get desired grid dimension
dim = nargin;

G.ni = size(x,2);
if dim>1, G.nj = size(y,2); else G.nj = 1; end
if dim>2, G.nk = size(z,2); else G.nk = 1; end

G.nn = G.ni * G.nj * G.nk;

if dim==3,
    G.ne = (G.ni-1) * (G.nj-1) * (G.nk-1);
elseif dim==2,
    G.ne = (G.ni-1) * (G.nj-1); 
else
    G.ne = G.ni-1;
end

G.dimension = dim;
if     dim==3, G.orientation = '3D';
elseif dim==2, G.orientation = '2D';
else           G.orientation = '1D';
end

if dim==3,
    [X,Y,Z] = ndgrid(x,y,z);  
    G.coord = zeros(G.nn,3);
    G.coord(:,1) = reshape(X,G.nn,1);
    G.coord(:,2) = reshape(Y,G.nn,1);
    G.coord(:,3) = reshape(Z,G.nn,1);
    
elseif dim==2,
    [X,Y] = ndgrid(x,y);
    G.coord = zeros(G.nn,2);
    G.coord(:,1) = reshape(X,G.nn,1);
    G.coord(:,2) = reshape(Y,G.nn,1);
    
else
    G.coord = zeros(G.nn,1);
    G.coord(:,1) = x;
    
end

if dim==3,     nn_e = 8;
elseif dim==2, nn_e = 4;
elseif dim==1, nn_e = 2;
end
G.elem_nodes = zeros(G.ne,nn_e);

ni = G.ni; nj = G.nj; nk = G.nk;
if dim==3
   for k=1:nk-1
       for j=1:nj-1
           for i=1:ni-1
               ie = i + (j-1)*(ni-1) + (k-1)*(ni-1)*(nj-1);
               G.elem_nodes(ie,1) = i + (j-1)*ni + (k-1)*ni*nj;
               G.elem_nodes(ie,2) = G.elem_nodes(ie,1) + 1;
               G.elem_nodes(ie,3) = G.elem_nodes(ie,2) + ni;
               G.elem_nodes(ie,4) = G.elem_nodes(ie,3) - 1;
               G.elem_nodes(ie,5) = G.elem_nodes(ie,1) + ni*nj;
               G.elem_nodes(ie,6) = G.elem_nodes(ie,2) + ni*nj;
               G.elem_nodes(ie,7) = G.elem_nodes(ie,3) + ni*nj;
               G.elem_nodes(ie,8) = G.elem_nodes(ie,4) + ni*nj;
           end
       end
   end
elseif dim==2
       for j=1:nj-1
           for i=1:ni-1
               ie = i + (j-1)*(ni-1);
               G.elem_nodes(ie,1) = i + (j-1)*ni;
               G.elem_nodes(ie,2) = G.elem_nodes(ie,1) + 1;
               G.elem_nodes(ie,3) = G.elem_nodes(ie,2) + ni;
               G.elem_nodes(ie,4) = G.elem_nodes(ie,3) - 1;
           end
       end
end

end