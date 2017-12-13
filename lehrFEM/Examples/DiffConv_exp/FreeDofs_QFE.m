function FreeDofs = FreeDofs_QFE(Mesh,BdFlags,varargin)
% FreeDofs_QFE free Dofs.
%
%   FREEDOFS is a M-by-1 matrix specifying the vertices and edges with no
%   prescribed Dirichlet boundary data.  
%
%   Example:
%
%   FreeDofs =FreeDofs_QFE(Mesh,BdFlags);
%
%   See also get_BdEdges.

%   Copyright 2005-2008 Patrick Meury, Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants

  nCoordinates = size(Mesh.Coordinates,1);
  nEdges = size(Mesh.Edges,1);
  
  tmp_1 = [];
  tmp_2 = [];
  U = zeros(nCoordinates+nEdges,1);
  for j = BdFlags
  
    % Extract Dirichlet nodes
  
    Loc = get_BdEdges(Mesh);
    DEdges = Loc(Mesh.BdFlags(Loc) == j);
    DNodes = unique([Mesh.Edges(DEdges,1); Mesh.Edges(DEdges,2)]);
                   
    % Collect Dirichlet nodes in temporary container
    
    tmp_1 = [tmp_1; DNodes];
    tmp_2 = [tmp_2; DEdges];
    
  end
  
  % Compute set of free dofs
  
  FreeDofs = [setdiff(1:nCoordinates,tmp_1) ...
              setdiff(1:nEdges,tmp_2) + nCoordinates];
  
return