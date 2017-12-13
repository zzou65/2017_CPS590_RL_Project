function FreeDofs = FreeDofs_LFE(Mesh,BdFlags,varargin)
% FreeDofs_LFE free DOFS 
%
%   FREEDOFS is a M-by-1 matrix specifying the vertices with no prescribed
%   Dirichlet boundary data.  
%
%   Example:
%
%   FreeDofs = FreeDofs_LFE(Mesh,BdFlags);
%
%   See also get_BdEdges.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants

  nCoordinates = size(Mesh.Coordinates,1);
  tmp = [];
 
  for j = BdFlags
  
    % Extract Dirichlet nodes
  
    Loc = get_BdEdges(Mesh);
    DEdges = Loc(Mesh.BdFlags(Loc) == j);
    DNodes = unique([Mesh.Edges(DEdges,1); Mesh.Edges(DEdges,2)]);
  
    % Collect Dirichlet nodes in temporary container
    
    tmp = [tmp; DNodes];
    
  end
  
  % Compute set of free dofs
  
  FreeDofs = setdiff(1:nCoordinates,tmp);
  
return