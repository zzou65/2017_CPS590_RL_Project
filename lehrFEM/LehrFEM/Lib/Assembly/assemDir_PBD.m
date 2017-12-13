function [U,FreeDofs] = assemDir_PBD(Mesh,BdFlags,FHandle,varargin)
% ASSEMDIR_PBD Dirichlet boundary conditions.
%
%   [U,FREEDOFS] = ASSEMDIR_PBD(MESH,BDFLAG,FHANDLE) incoporates the
%   Dirichlet boundary conditions with the data given by FHANDLE into the
%   finite element solution U. The boundary condition is only enforced at
%   the vertices of the edges whose BdFlag is equal to the integer BDFLAG.
%
%   [U,FREEDOFS] = ASSEMDIR_PBD(MESH,BDFLAG,FHANDLE,FPARAM) also
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
%   [U,FreeDofs] = assemDir_PBD(Mesh,BdFlags,FHandle);
%
%   See also get_BdEdges.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants

  Rot = [0 -1; 1 0];
  nCoordinates = size(Mesh.Coordinates,1);
  nEdges = size(Mesh.Edges,1);
  
  tmp_1 = [];
  tmp_2 = [];
  U = zeros(nCoordinates+nEdges,1);
  for j1 = BdFlags
  
    % Extract boundary edges, vertices and elements
  
    Loc = get_BdEdges(Mesh);
    DEdges = Loc(Mesh.BdFlags(Loc) == j1);
    DNodes = unique([Mesh.Edges(DEdges,1); Mesh.Edges(DEdges,2)]);
    [rows,cols,vals] = find(Mesh.Edge2Elem(DEdges,:));
    vals = sortrows([rows cols vals],1);
    DElem = vals(:,3);
    
    % Compute mid points of curved edges
   
    normal = Mesh.Coordinates(Mesh.Elements(DElem,3),:) - ...
             Mesh.Coordinates(Mesh.Elements(DElem,2),:);         
    normal = (normal*Rot)./(sqrt(sum(normal.^2,2))*ones(1,2));
    XMid = (Mesh.Coordinates(Mesh.Elements(DElem,2),:) + ...
            Mesh.Coordinates(Mesh.Elements(DElem,3),:))/2 + ...
           (Mesh.Delta(DEdges)*ones(1,2)).*normal;
           
    % Incoporate Dirichlet boundary conditions
    
    U(DNodes) = FHandle(Mesh.Coordinates(DNodes,:),j,varargin{:});  
    U(DEdges+nCoordinates) = FHandle(XMid,j,varargin{:});
    
    % Collect Dirichlet nodes in temporary container
    
    tmp_1 = [tmp_1; DNodes];
    tmp_2 = [tmp_2; DEdges]; 
    
  end
  
  % Compute set of free dofs
 
  FreeDofs = [setdiff(1:nCoordinates,tmp_1) ...
              setdiff(1:nEdges,tmp_2) + nCoordinates];
          
return