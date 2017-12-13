function varargout = convert_LFE2_W1F(Mesh)
% CONVERT_LFE2_W1F transformation matrix that convert Nodal FE solution to 
% Edge FE solution.
%
%   P = CONVERT_LFE2_W1F(MESH) gives us the transformation matrix P.
%  
%   P is a P-by-2M matrix that can do the following convertion:
%   U_W1F = P * U_LFE2
%
%   [I,J,P] = CONVERT_LFE2_W1F(MESH) returns the matrix in an array 
%   representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 or N-by-4 matrix specifying the elements of the
%                 mesh.
%    EDGES        P-by-2 matrix specifying all edges of the mesh.
%
%   Example:
%
%   P = convert_LFE2_W1F(Mesh);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    nEdges = size(Mesh.Edges,1);
    nCoordinates = size(Mesh.Coordinates,1);
    
    % Preallocate memory
    
    I = zeros(4*nEdges,1); 
    J = zeros(4*nEdges,1); 
    P = zeros(4*nEdges,1); 
    
    % Loop all the edges
    
    loc = 1:4;
    for i = 1:nEdges
       
       vidx = Mesh.Edges(i,:);
       
       % Compute tangential vector
       
       t = 1/2*(Mesh.Coordinates(vidx(2),:)-Mesh.Coordinates(vidx(1),:));
       
       % Add contributions to the matrix
       
       tI = [i;i;i;i];
       tJ = [vidx(1);vidx(2);vidx(1)+nCoordinates;vidx(2)+nCoordinates];
       tP = [t(1);t(1);t(2);t(2)];
       
       I(loc) = tI;
       J(loc) = tJ;
       P(loc) = tP;
       loc = loc + 4;
       
    end
    
    % Assign output arguments
    
    if(nargout > 1)
        varargout{1} = I;
        varargout{2} = J;
        varargout{3} = P;
    else
        varargout{1} = sparse(I,J,P);
    end

return