function [U,FreeDofs] = assemDir_Stokes_TH(Mesh,BdFlags,FHandle,varargin)
% ASSEMDIR_STOKES_TH Dirichlet boundary conditions for Taylor-Hood elements.
%
%   [U,FREEDOFS] = ASSEMDIR_STOKES_TH(MESH,BDFLAG,FHANDLE) incoporates
%   the Dirichlet boundary conditions with the data given by FHANDLE into
%   the finite element solution U. The boundary condition is only enforced
%   at the vertices of the edges whose BdFlag is equal to the integer
%   BDFLAG.
%
%   [U,FREEDOFS] = ASSEMDIR_STOKES_TH(MESH,BDFLAG,FHANDLE,FPARAM) also
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
%   [U,FreeDofs] = assemDir_Stokes_TH(Mesh,BdFlags,FHandle);
%
%   See also get_BdEdges.

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants

  nCoordinates = size(Mesh.Coordinates,1);
  nEdges = size(Mesh.Edges,1);
  nElements = size(Mesh.Elements,1);
  
  tmp_1 = [];
  tmp_2 = [];
  U = zeros(3*nCoordinates+2*nEdges+1,1);
  for j = BdFlags
  
    % Extract Dirichlet nodes
  
    Loc = get_BdEdges(Mesh);
    DEdges = Loc(Mesh.BdFlags(Loc) == j);
    DNodes = unique([Mesh.Edges(DEdges,1); Mesh.Edges(DEdges,2)]);
    
    % Compute midpoints of all edges
    
    MidPoints = 1/2*(Mesh.Coordinates(Mesh.Edges(DEdges,1),:) ...
                   + Mesh.Coordinates(Mesh.Edges(DEdges,2),:));
    
    % Compute Dirichlet boundary conditions
  
    FVal = FHandle(Mesh.Coordinates(DNodes,:),j,varargin{:});
    U(DNodes) = FVal(:,1);
    U(DNodes+nCoordinates+nEdges) = FVal(:,2);
    
    FVal = FHandle(MidPoints,j,varargin{:});
    U(DEdges+nCoordinates) = FVal(:,1);
    U(DEdges+2*nCoordinates+nEdges) = FVal(:,2);               
    
    % Collect Dirichlet nodes in temporary container
    
    tmp_1 = [tmp_1; DNodes];
    tmp_2 = [tmp_2; DEdges];
    
  end
  
  % Compute free dofs
  
  z1 = setdiff(1:nCoordinates,unique(tmp_1));
  z2 = setdiff(1:nEdges,unique(tmp_2));
  FreeDofs = [z1 ...
              z2+nCoordinates ...
              z1+nCoordinates+nEdges ...
              z2+2*nCoordinates+nEdges ...
              2*(nCoordinates+nEdges)+(1:nCoordinates) ...
              3*nCoordinates+2*nEdges+1];
    
return