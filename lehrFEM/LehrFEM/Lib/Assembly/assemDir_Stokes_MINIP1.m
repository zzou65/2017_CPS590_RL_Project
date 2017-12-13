function [U,FreeDofs] = assemDir_Stokes_MINIP1(Mesh,BdFlags,FHandle,varargin)
% ASSEMDIR_STOKES_MINIP1 Dirichlet boundary conditions.
%
%   [U,FREEDOFS] = ASSEMDIR_STOKES_MINIP1(MESH,BDFLAG,FHANDLE) incoporates
%   the Dirichlet boundary conditions with the data given by FHANDLE into
%   the finite element solution U. The boundary condition is only enforced
%   at the vertices of the edges whose BdFlag is equal to the integer
%   BDFLAG.
%
%   [U,FREEDOFS] = ASSEMDIR_STOKES_MINIP1(MESH,BDFLAG,FHANDLE,FPARAM) also
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
%   [U,FreeDofs] = assemDir_Stokes_MINIP1(Mesh,BdFlags,FHandle);
%
%   See also get_BdEdges.

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants

  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  
  tmp = [];
  U = zeros(3*nCoordinates+2*nElements+1,1);
  for j = BdFlags
  
    % Extract Dirichlet nodes
  
    Loc = get_BdEdges(Mesh);
    DEdges = Loc(Mesh.BdFlags(Loc) == j);
    DNodes = unique([Mesh.Edges(DEdges,1); Mesh.Edges(DEdges,2)]);
  
    % Compute Dirichlet boundary conditions
  
    FVal = FHandle(Mesh.Coordinates(DNodes,:),j,varargin{:});
    U(DNodes) = FVal(:,1);
    U(DNodes+nCoordinates+nElements) = FVal(:,2);
    
    % Collect Dirichlet nodes for each component in temporary containers
    
    tmp = [tmp; DNodes];
        
  end
  
  % Compute free dofs
 
  z = setdiff(1:nCoordinates,unique(tmp));
  FreeDofs = [z ...
              (1:nElements)+nCoordinates ...
              z+nCoordinates+nElements ...
              (1:nElements)+2*nCoordinates+nElements ...
              (1:nCoordinates)+2*(nCoordinates+nElements) ...
              3*nCoordinates+2*nElements+1];
    
return