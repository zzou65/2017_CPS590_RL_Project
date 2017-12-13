function [U,FreeDofs] = assemDir_Stokes_CRP0(Mesh,BdFlags,FHandle,varargin)
% ASSEMDIR_STOKES_CRP0 Dirichlet boundary conditions.
%
%   [U,FREEDOFS] = ASSEMDIR_STOKES_CRP0(MESH,BDFLAG,FHANDLE) incoporates
%   the Dirichlet boundary conditions with the data given by FHANDLE into
%   the finite element solution U. The boundary condition is only enforced
%   at the vertices of the edges whose BdFlag is equal to the integer
%   BDFLAG.
%
%   [U,FREEDOFS] = ASSEMDIR_STOKES_CRP0(MESH,BDFLAG,FHANDLE,FPARAM) also
%   handles the variable length argument list FPARAM to the boundary data
%   function FHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    EDGES        P-by-2 matrix specifying the edges of the mesh.
%    BDFLAGS      P-by-1 matrix specifying the boundary type of each 
%                 boundary edge in the mesh.
%
%   FREEDOFS is a M-by-1 matrix specifying the dofs with no prescribed
%   Dirichlet boundary data.  
%
%   Example:
%
%   [U,FreeDofs] = assemDir_Stokes_CRP0(Mesh,BdFlags,FHandle);
%
%   See also get_BdEdges.

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants

  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);
  
  tmp = [];
  U = zeros(2*nEdges+nElements+1,1);
  for j = BdFlags
  
    % Extract Dirichlet nodes
  
    Loc = get_BdEdges(Mesh);
    Loc = Loc(Mesh.BdFlags(Loc) == j);

    % Compute Dirichlet boundary conditions

    MidPoints = 1/2*(Mesh.Coordinates(Mesh.Edges(Loc,1),:) + ...
                     Mesh.Coordinates(Mesh.Edges(Loc,2),:));
    FVal = FHandle(MidPoints,j,varargin{:});
    U(Loc) = FVal(:,1);
    U(Loc+nEdges) = FVal(:,2);
    
    % Collect Dirichlet edges in temporary container
    
    tmp = [tmp; Loc];
    
  end
  
  % Eliminate mulitple vertices
  
  z = setdiff(1:nEdges,unique(tmp));
  FreeDofs = [z ...
              z+nEdges ...
              2*nEdges+(1:nElements) ...
              2*nEdges+nElements+1];
  
return