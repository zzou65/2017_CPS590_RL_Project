function [U,FreeDofs] = assemDir_LFE(Mesh,BdFlags,FHandle,varargin)
% ASSEMDIR_LFE Dirichlet boundary conditions.
%
%   [U,FREEDOFS] = ASSEMDIR_LFE(MESH,BDFLAG,FHANDLE) incoporates the
%   Dirichlet boundary conditions with the data given by FHANDLE into the
%   finite element solution U. The boundary condition is only enforced at
%   the vertices of the edges whose BdFlag is equal to the integer BDFLAG.
%
%   [U,FREEDOFS] = ASSEMDIR_LFE(MESH,BDFLAG,FHANDLE,FPARAM) also
%   handles the variable length argument list FPARAM to the boundary data
%   function FHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    EDGES        P-by-2 matrix specifying the edges of the mesh.
%    BDFLAGS      P-by-1 matrix specifying the boundary type of each
%                 boundary edge in the mesh.
%
%   FREEDOFS is a M-by-1 matrix specifying the vertices with no prescribed
%   Dirichlet boundary data.  
%
%   Example:
%
%   [U,FreeDofs] = assemDir_LFE(Mesh,BdFlags,FHandle);
%
%   See also get_BdEdges.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants

  nCoordinates = size(Mesh.Coordinates,1);
  tmp = [];
  U = zeros(nCoordinates,1);
 
  for j = BdFlags
  
    % Extract Dirichlet nodes
  
    Loc = get_BdEdges(Mesh);
    DEdges = Loc(Mesh.BdFlags(Loc) == j);
    DNodes = unique([Mesh.Edges(DEdges,1); Mesh.Edges(DEdges,2)]);
  
    % Compute Dirichlet boundary conditions
  
    U(DNodes) = FHandle(Mesh.Coordinates(DNodes,:),j,varargin{:});  
  
    % Collect Dirichlet nodes in temporary container
    
    tmp = [tmp; DNodes];
    
  end
  
  % Compute set of free dofs
  
  FreeDofs = setdiff(1:nCoordinates,tmp);
  
return