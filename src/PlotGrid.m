function PlotGrid(G, varargin) 
% Plots the mesh or scalar data on external (visible) faces of a 2D/3D grid 
%
% INPUTS:
% G                 - Cartesian grid structure as created by calling 
%                     cartGrid routine   
% varargin (optional)
%                   - optional property-value pair arguments submitted to 
%                     MATLAB patch command. Type 'help patch' in the 
%                     command window to see the list of supported options
%
% Author: M.A. Sbai, Ph.D.
%         BRGM (French Geological Survey) 
%         D3E  (Direction Eau, Environnement, Echotechnologies)
% 

ni = G.ni;
nj = G.nj; 
nk = G.nk;

blanked = false;
if ~isempty(varargin) && isa(varargin{1}, 'double') 
    blanked = true;
    blank  = varargin{1};
    assert(size(blank,1)==3 && size(blank,2)==2, ...
        'Size of blank zone must be 3x2');
    varargin(1) = [];
end

% implementation for 2D case  
if strcmp(G.orientation,'2D'), 
    
    num = zeros(nj,ni);
    for j=1:nj
        for i=1:ni
            num(j,i) = i + (j-1)*ni;
        end
    end

    coord = zeros(4,(ni-1)*(nj-1),2);
    c = 1;
    for j=1:nj-1
        for i=1:ni-1
            coord(1,c,:) = G.coord(num(j,i),:);
            coord(2,c,:) = G.coord(num(j,i+1),:);
            coord(3,c,:) = G.coord(num(j+1,i+1),:);
            coord(4,c,:) = G.coord(num(j+1,i),:);
            c = c + 1; 
        end
    end
    
    patch('XData',coord(:,:,1), 'YData',coord(:,:,2), varargin{:});
    
    return
end 

% --- Plot the six faces of the 3D mesh ---
% 

i1b = false; i2b = false; 
j1b = false; j2b = false; 
k1b = false; k2b = false;
if blanked 
    cartDims = [ni nj nk];
    i1b = isBlanked(cartDims, blank, true,  'I');
    i2b = isBlanked(cartDims, blank, false, 'I');
    j1b = isBlanked(cartDims, blank, true,  'J');
    j2b = isBlanked(cartDims, blank, false, 'J');
    k1b = isBlanked(cartDims, blank, true,  'K');
    k2b = isBlanked(cartDims, blank, false, 'K');
end

% (1) Left  face @ I=1
if ~i1b 
    coord = extractSliceNodes(G, 1, 'I');  
else
    coord = extractSliceNodes(G, 1, 'I', blank(2:3,:));
end
patch(coord(:,:,1), coord(:,:,2), coord(:,:,3), varargin{:});

% (2) right face @ I=NI
if ~i2b
    coord = extractSliceNodes(G, ni, 'I');
else
    coord = extractSliceNodes(G, ni, 'I', blank(2:3,:));
end
patch(coord(:,:,1), coord(:,:,2), coord(:,:,3), varargin{:});

% (3) front face @ J=1
if ~j1b
    coord = extractSliceNodes(G, 1, 'J');
else
    coord = extractSliceNodes(G, 1, 'J', [blank(1,:);blank(3,:)]);
end
patch(coord(:,:,1), coord(:,:,2), coord(:,:,3), varargin{:});

% (4) back  face @ J=NJ 
if ~j2b
    coord = extractSliceNodes(G, nj, 'J'); 
else
    coord = extractSliceNodes(G, nj, 'J', [blank(1,:);blank(3,:)]); 
end
patch(coord(:,:,1), coord(:,:,2), coord(:,:,3), varargin{:});

% (5) top  face @ K=1
if ~k1b
    coord = extractSliceNodes(G, 1, 'K');  
else
    coord = extractSliceNodes(G, 1, 'K', blank(1:2,:));  
end
patch(coord(:,:,1), coord(:,:,2), coord(:,:,3), varargin{:});

% (6) bottom face @ K=NK
if ~k2b
    coord = extractSliceNodes(G, nk, 'K'); 
else
    coord = extractSliceNodes(G, nk, 'K', blank(1:2,:)); 
end
patch(coord(:,:,1), coord(:,:,2), coord(:,:,3), varargin{:});
 
% --- Plot possible external faces after grid splitting along given slices 
% 
if isfield(G, 'split_nodes')
     
    num_ss = length(G.split_nodes);         % number of split slices
    
    for is=1:num_ss
        k  = cell2mat(G.split_nodes{is}(1));    % slice index  
        ij = cell2mat(G.split_nodes{is}(2));    % split nodes for k-slice 
        
        up = true;
        up_elems = neighborElements(G, ij, k, up);
        up_elems = unique(up_elems,'stable');
        coord = extractElemFaceNodes(G, up_elems, ~up);
        
        patch(coord(:,:,1), coord(:,:,2), coord(:,:,3), varargin{:});
        
        up = false;
        down_elems = neighborElements(G, ij, k, up);
        down_elems = unique(down_elems,'stable');
        coord = extractElemFaceNodes(G, down_elems, ~up);
        
        patch(coord(:,:,1), coord(:,:,2), coord(:,:,3), varargin{:});
    end
    
end 

% determine the number of other I, J, K patches to plot 
%%nf = 6 - sum([i1b i2b j1b j2b k1b k2b]);

if blanked 
    
    if blank(1,1)>1  % new plot 
        coord = extractSliceNodes(G, blank(1,1), 'I', blank(2:3,:), true);
        patch(coord(:,:,1), coord(:,:,2), coord(:,:,3), varargin{:});
    end

    if blank(1,2)<ni
        coord = extractSliceNodes(G, blank(1,2), 'I', blank(2:3,:), true);
        patch(coord(:,:,1), coord(:,:,2), coord(:,:,3), varargin{:});
    end

    if blank(2,1)>1 
        coord = extractSliceNodes(G, blank(2,1), 'J', [blank(1,:);blank(3,:)], true);
        patch(coord(:,:,1), coord(:,:,2), coord(:,:,3), varargin{:});
    end

    if blank(2,2)<nj
        coord = extractSliceNodes(G, blank(2,2), 'J', [blank(1,:);blank(3,:)], true);
        patch(coord(:,:,1), coord(:,:,2), coord(:,:,3), varargin{:});
    end

    if blank(3,1)>1 
        coord = extractSliceNodes(G, blank(3,1), 'K', blank(1:2,:), true);
        patch(coord(:,:,1), coord(:,:,2), coord(:,:,3), varargin{:});
    end

    if blank(3,2)<nk 
        coord = extractSliceNodes(G, blank(3,2), 'K', blank(1:2,:), true);
        patch(coord(:,:,1), coord(:,:,2), coord(:,:,3), varargin{:});
    end

end

view(3), axis equal 
end 


function b = isBlanked(cartDims, blank, first, direction) 

assert(length(cartDims)==3, 'Size of cartDims must be 1x3');
assert(first==true || first==false, 'First argument must be true or false.');
assert(strcmp(direction,'I') || strcmp(direction,'J') || ... 
       strcmp(direction,'K'), 'direction must be either I, J, or K.');

b = false;
if strcmp(direction,'I') && first==true
    b = (blank(1,1)==1);
end

if strcmp(direction,'I') && first==false
   b = (blank(1,2)==cartDims(1)); 
end

if strcmp(direction,'J') && first==true
    b = (blank(2,1)==1);
end

if strcmp(direction,'J') && first==false
   b = (blank(2,2)==cartDims(2)); 
end

if strcmp(direction,'K') && first==true
    b = (blank(3,1)==1);
end

if strcmp(direction,'K') && first==false
   b = (blank(3,2)==cartDims(3)); 
end
   
end


function elems = neighborElements(G, ij, k, up) 
% Returns indices of elements holding the nodes given by indices array ij 
% and in slice k from top (up = true) or down (up = false) directions. 
% 

assert(size(ij,2)==2);
assert(k>=1 || k<=G.nk);
assert(isa(up,'logical'));

ni = G.ni; nj = G.nj; nk = G.nk;
num_nodes = size(ij,1);
elems = [];

if ~up 
    if k==nk, return; end
    for n=1:num_nodes
        i = ij(n,1);
        j = ij(n,2);
        if i>1  && j>1,  elems = [elems element(G,[i-1 j-1 k])]; end
        if i<ni && j>1,  elems = [elems element(G,[i   j-1 k])]; end
        if i>1  && j<nj, elems = [elems element(G,[i-1 j   k])]; end
        if i<ni && j<nj, elems = [elems element(G,[i   j   k])]; end
    end
else
    if k==1, return; end
    for n=1:num_nodes
        i = ij(n,1);
        j = ij(n,2);
        if i>1  && j>1,  elems = [elems element(G,[i-1 j-1 k-1])]; end
        if i<ni && j>1,  elems = [elems element(G,[i   j-1 k-1])]; end
        if i>1  && j<nj, elems = [elems element(G,[i-1 j   k-1])]; end
        if i<ni && j<nj, elems = [elems element(G,[i   j   k-1])]; end
    end
end

end