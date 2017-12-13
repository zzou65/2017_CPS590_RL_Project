function varargout = MAT_G2E(Mesh)
% MAT_G2E transformation matrix from LFE grad space to edge element space
%
%   GE = MAT_G2E(MESH) gives us the transformation matrix GE.
%  
%   GE is a P-by-M matrix that can do the following convertion:
%   grad:{LFE(0 on BD)} -> Edge elements
%
%   [I,J,GE] = MAT_G2E(MESH) returns the matrix in an array representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 or N-by-4 matrix specifying the elements of the
%                 mesh.
%    EDGES        P-by-2 matrix specifying all edges of the mesh.
%
%   Example:
%
%   GE = Mat_G2E(Mesh);

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    nEdges = size(Mesh.Edges,1);
    nCoordinates = size(Mesh.Coordinates,1);
    
    % Preallocate memory
    
    I = zeros(2*nEdges,1); 
    J = zeros(2*nEdges,1); 
    G = zeros(2*nEdges,1); 
    
    % Loop all the edges
    
    loc = 1:2;
    for i = 1:nEdges
       
       vidx = Mesh.Edges(i,:);

       % Add contributions to the matrix
       
       tI = [i;i];
       tJ = [vidx(1);vidx(2)];
       tP = [-1;1];
       
       I(loc) = tI;
       J(loc) = tJ;
       G(loc) = tP;
       loc = loc + 2;
       
    end
    
    % Assign output arguments
    
    if(nargout > 1)
        varargout{1} = I;
        varargout{2} = J;
        varargout{3} = G;
    else
        varargout{1} = sparse(I,J,G);
    end

return