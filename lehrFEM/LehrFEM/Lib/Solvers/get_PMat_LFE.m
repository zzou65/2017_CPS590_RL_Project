function P = get_PMat_LFE(CMesh,FMesh)
% GET_PMAT_LFE Compute prolongation matrix.
%
%   P = GET_PMAT_LFE(MESH) computes the matrix corresponding to the
%   prolongation operator from a coarse mesh CMESH to a fine mesh FMESH
%   obtained through regular red refinements using REFINE_REG or
%   REFINE_REG_JIGGLE.
%
%   The structs CMESH and FMESH should contain at least the following
%   fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh. 
%    EDGES       P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE   M-by-M sparse matrix which specifies wheter the two vertices
%                i and j are connected by an edge with number VERT2EDGE(i,j).
%   
%   Example:
%
%   P = get_PMat_LFE(CMesh,FMesh);
%
%   See also add_MLevel.

%   Copyright 2005-2005 Patrick Meury & Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nCCoordinates = size(CMesh.Coordinates,1);
  nFCoordinates = size(FMesh.Coordinates,1);
  
  % Add multilevel information to the fine mesh
  
  FMesh = add_MLevel(FMesh);
  
  % Preallocate memory
  
  nEntries = nCCoordinates + 2*(nFCoordinates-nCCoordinates);
  I = zeros(nEntries,1);
  J = zeros(nEntries,1);
  P = zeros(nEntries,1);
  
  % Build prolongation matrix
 
  nnz = 0;
  for i = 1:nFCoordinates
    if(FMesh.Father_Vert(i,1) > 0)
      I(nnz+1) = i;
      J(nnz+1) = FMesh.Father_Vert(i,1);
      P(nnz+1) = 1;
      nnz = nnz+1;
    else  
      I(nnz+1) = i;
      I(nnz+2) = i;
      J(nnz+1) = CMesh.Edges(FMesh.Father_Vert(i,2),1);
      J(nnz+2) = CMesh.Edges(FMesh.Father_Vert(i,2),2);
      P(nnz+1) = norm(FMesh.Coordinates(i,:)-CMesh.Coordinates(J(nnz+2),:))...
          /norm(CMesh.Coordinates(J(nnz+1),:)-CMesh.Coordinates(J(nnz+2),:));
      P(nnz+2) = 1-P(nnz+1);
      nnz = nnz+2;
    end
  end
  
  % Convert to sparse matrix
  
  P = sparse(I,J,P);
  
return