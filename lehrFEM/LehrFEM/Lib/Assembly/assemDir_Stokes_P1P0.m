function [U,FreeDofs] = assemDir_Stokes_P1P0(Mesh,BdFlags,FHandle,varargin)
% ASSEMDIR_STOKES_P1P0 Dirichlet boundary conditions for P1P0 elements.
%
%   [U,FREEDOFS] = ASSEMDIR_STOKES_P1P0(MESH,BDFLAG,FHANDLE) incoporates
%   the Dirichlet boundary conditions with the data given by FHANDLE into
%   the finite element solution U. The boundary condition is only enforced
%   at the vertices of the edges whose BdFlag is equal to the integer
%   BDFLAG.
%
%   [U,FREEDOFS] = ASSEMDIR_STOKES_P1P0(MESH,BDFLAG,FHANDLE,FPARAM) also
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
%   [U,FreeDofs] = assemDir_Stokes_P1P0(Mesh,BdFlags,FHandle);
%
%   See also get_BdEdges.

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants

  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  
  tmp = [];
  U = zeros(2*nCoordinates+nElements+1,1);
  for j = BdFlags
  
    % Extract Dirichlet nodes
  
    Loc = get_BdEdges(Mesh);
    DEdges = Loc(Mesh.BdFlags(Loc) == j);
    DNodes = unique([Mesh.Edges(DEdges,1); Mesh.Edges(DEdges,2)]);
    
    % Compute Dirichlet boundary conditions
  
    FVal = FHandle(Mesh.Coordinates(DNodes,:),j,varargin{:});
    U(DNodes) = FVal(:,1);
    U(DNodes+nCoordinates) = FVal(:,2);            
    
    % Collect Dirichlet nodes in temporary container
    
    tmp = [tmp; DNodes];
    
  end
  
  % Compute free dofs
  
  z1 = setdiff(1:nCoordinates,unique(tmp));
  FreeDofs = [z1 ...
              z1+nCoordinates ...
              2*nCoordinates+(1:nElements) ...
              2*nCoordinates+nElements+1];
    
return