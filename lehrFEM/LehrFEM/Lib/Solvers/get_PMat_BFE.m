function P = get_PMat_BFE(CMesh,FMesh)
%GET_PMAT_BFE Compute prolongation matrix for BFE
%
%   P = GET_PMAT_LFE(MESH) computes the matrix corresponding to the
%   prolongation operator from a coarse mesh CMESH to a fine mesh FMESH
%   obtained through regular red refinements using REFINE_REG or
%   REFINE_REG_JIGGLE.
%
%   The structs CMESH and FMESH should contain at least the following
%   fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-4 matrix specifying the elements of the mesh. 
%    EDGES       P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE   M-by-M sparse matrix which specifies wheter the two vertices
%                i and j are connected by an edge with number VERT2EDGE(i,j).
%   
%   Example:
%
%   P = get_PMat_LFE(CMesh,FMesh);
%
%   See also add_MLevel, get_PMat_LFE.

%   Copyright 2006-2006 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nCCoordinates = size(CMesh.Coordinates,1);
  nCEdges = size(CMesh.Edges,1);
  nCElements = size(CMesh.Elements,1);
  nFCoordinates = size(FMesh.Coordinates,1);
  
  % Add multilevel information to the fine mesh
  
  FMesh = add_MLevel(FMesh);
  
  % Preallocate memory
  
  nEntries = nCCoordinates + 2*nCEdges + 4*nCElements;
  I = zeros(nEntries,1);
  J = zeros(nEntries,1);
  P = zeros(nEntries,1);
  
  nnz = 0;
  for i = 1:nFCoordinates
    % i is vertex of CMesh
    if(FMesh.Father_Vert(i,1) > 0)
      I(nnz+1) = i;
      J(nnz+1) = FMesh.Father_Vert(i,1);
      P(nnz+1) = 1;
      nnz = nnz+1;
    % i is on edge of CMesh
    elseif(FMesh.Father_Vert(i,2) > 0)
      I(nnz+1) = i;
      I(nnz+2) = i;
      J(nnz+1) = CMesh.Edges(FMesh.Father_Vert(i,2),1);
      J(nnz+2) = CMesh.Edges(FMesh.Father_Vert(i,2),2);
      P(nnz+1) = norm(FMesh.Coordinates(i,:)-CMesh.Coordinates(J(nnz+2),:))...
          /norm(CMesh.Coordinates(J(nnz+1),:)-CMesh.Coordinates(J(nnz+2),:));
      P(nnz+2) = 1-P(nnz+1);
      nnz = nnz+2;
    % i is inside element of CMesh
    else
      I(nnz+1) = i;
      I(nnz+2) = i;
      I(nnz+3) = i;
      I(nnz+4) = i;
      J(nnz+1) = CMesh.Elements(FMesh.Father_Vert(i,3),1);
      J(nnz+2) = CMesh.Elements(FMesh.Father_Vert(i,3),2);
      J(nnz+3) = CMesh.Elements(FMesh.Father_Vert(i,3),3);
      J(nnz+4) = CMesh.Elements(FMesh.Father_Vert(i,3),4);
        
      B = [CMesh.Coordinates(J(nnz+2),:)-CMesh.Coordinates(J(nnz+1),:);...
        CMesh.Coordinates(J(nnz+4),:)-CMesh.Coordinates(J(nnz+1),:)];
      C = (FMesh.Coordinates(i,:)-CMesh.Coordinates(J(nnz+1),:))/B;
      P(nnz+1) = (1-C(1))*(1-C(2));
      
      B = [CMesh.Coordinates(J(nnz+3),:)-CMesh.Coordinates(J(nnz+2),:);...
        CMesh.Coordinates(J(nnz+1),:)-CMesh.Coordinates(J(nnz+2),:)];
      C = (FMesh.Coordinates(i,:)-CMesh.Coordinates(J(nnz+2),:))/B;
      P(nnz+2) = (1-C(1))*(1-C(2));
      
      B = [CMesh.Coordinates(J(nnz+4),:)-CMesh.Coordinates(J(nnz+3),:);...
        CMesh.Coordinates(J(nnz+2),:)-CMesh.Coordinates(J(nnz+3),:)];
      C = (FMesh.Coordinates(i,:)-CMesh.Coordinates(J(nnz+3),:))/B;
      P(nnz+3) = (1-C(1))*(1-C(2));
      
      B = [CMesh.Coordinates(J(nnz+1),:)-CMesh.Coordinates(J(nnz+4),:);...
        CMesh.Coordinates(J(nnz+3),:)-CMesh.Coordinates(J(nnz+4),:)];
      C = (FMesh.Coordinates(i,:)-CMesh.Coordinates(J(nnz+4),:))/B;
      P(nnz+4) = (1-C(1))*(1-C(2));
      
      nnz = nnz+4;
    end
  end
    
  % Convert to sparse matrix
  
  P = sparse(I,J,P);
  
return
  