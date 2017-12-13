function Mesh = add_MLevel(Mesh)
% ADD_MLEVEL Add multilevel information.
%
%   MESH = ADD_MLEVEL(MESH) adds mulitlevel information to the struct MESH.
%   It is assumed that the struct MESH was generated from a coarser grid by
%   uniform red refinements.
%
%   The structs MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the mesh. 
%    EDGES       P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE   M-by-M sparse matrix which specifies wheter the two vertices
%                i and j are connected by an edge with number VERT2EDGE(i,j).
%
%   Example:
%
%   Mesh = add_MLevel(Mesh);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  [nElements,nVert] = size(Mesh.Elements);
  nElements = nElements/4;
  if(nVert == 3)
    nEdges = (size(Mesh.Edges,1)-3*nElements)/2;
    nCoordinates = size(Mesh.Coordinates,1)-nEdges;   
      
    % Preallocate memory
  
    Mesh.Father_Elem = zeros(size(Mesh.Elements,1),1);
    Mesh.Father_Edge = zeros(size(Mesh.Edges,1),2);
    Mesh.Father_Vert = zeros(size(Mesh.Coordinates,1),3);
  
    for i = 1:nElements
      i1 = Mesh.Elements(4*(i-1)+1,1);
      i2 = Mesh.Elements(4*(i-1)+2,2);
      i3 = Mesh.Elements(4*(i-1)+3,3);
      j1 = Mesh.Elements(4*(i-1)+4,1);
      j2 = Mesh.Elements(4*(i-1)+4,2);
      j3 = Mesh.Elements(4*(i-1)+4,3);
      
      % Add father relationship for elements  
      
      Mesh.Father_Elem(4*(i-1)+1) = i;
      Mesh.Father_Elem(4*(i-1)+2) = i;
      Mesh.Father_Elem(4*(i-1)+3) = i;
      Mesh.Father_Elem(4*(i-1)+4) = i;
    
      % Add father relationship for vertices
    
      Mesh.Father_Vert(i1,1) = i1;
      Mesh.Father_Vert(i2,1) = i2;
      Mesh.Father_Vert(i3,1) = i3;
      Mesh.Father_Vert(j1,2) = j1-nCoordinates;
      Mesh.Father_Vert(j2,2) = j2-nCoordinates;
      Mesh.Father_Vert(j3,2) = j3-nCoordinates;
    
      % Add father relationship to edges
    
      Mesh.Father_Edge(Mesh.Vert2Edge(i2,j1),1) = j1-nCoordinates;
      Mesh.Father_Edge(Mesh.Vert2Edge(j1,i3),1) = j1-nCoordinates;
      Mesh.Father_Edge(Mesh.Vert2Edge(i3,j2),1) = j2-nCoordinates;
      Mesh.Father_Edge(Mesh.Vert2Edge(j2,i1),1) = j2-nCoordinates;
      Mesh.Father_Edge(Mesh.Vert2Edge(i1,j3),1) = j3-nCoordinates;
      Mesh.Father_Edge(Mesh.Vert2Edge(j3,i2),1) = j3-nCoordinates;
      Mesh.Father_Edge(Mesh.Vert2Edge(j1,j2),2) = i;
      Mesh.Father_Edge(Mesh.Vert2Edge(j2,j3),2) = i;
      Mesh.Father_Edge(Mesh.Vert2Edge(j3,j1),2) = i;      
    end

  else
    nEdges = (size(Mesh.Edges,1)-4*nElements)/2;
    nCoordinates = size(Mesh.Coordinates,1)-nEdges-nElements;   
    
    % Preallocate memory
  
    Mesh.Father_Elem = zeros(size(Mesh.Elements,1),1);
    Mesh.Father_Edge = zeros(size(Mesh.Edges,1),2);
    Mesh.Father_Vert = zeros(size(Mesh.Coordinates,1),3);
  
    for i = 1:nElements
      i1 = Mesh.Elements(4*(i-1)+1,1);
      i2 = Mesh.Elements(4*(i-1)+2,2);
      i3 = Mesh.Elements(4*(i-1)+4,3);
      i4 = Mesh.Elements(4*(i-1)+3,4);
      j1 = Mesh.Elements(4*(i-1)+1,2);
      j2 = Mesh.Elements(4*(i-1)+2,3);
      j3 = Mesh.Elements(4*(i-1)+3,3);
      j4 = Mesh.Elements(4*(i-1)+3,1);
      jc = Mesh.Elements(4*(i-1)+4,1);
      
      % Add father relationship for elements  
      
      Mesh.Father_Elem(4*(i-1)+1) = i;
      Mesh.Father_Elem(4*(i-1)+2) = i;
      Mesh.Father_Elem(4*(i-1)+3) = i;
      Mesh.Father_Elem(4*(i-1)+4) = i;
    
      % Add father relationship for vertices
    
      Mesh.Father_Vert(i1,1) = i1;
      Mesh.Father_Vert(i2,1) = i2;
      Mesh.Father_Vert(i3,1) = i3;
      Mesh.Father_Vert(i4,1) = i4;
      Mesh.Father_Vert(j1,2) = j1-nCoordinates;
      Mesh.Father_Vert(j2,2) = j2-nCoordinates;
      Mesh.Father_Vert(j3,2) = j3-nCoordinates;
      Mesh.Father_Vert(j4,2) = j4-nCoordinates;
      Mesh.Father_Vert(jc,3) = i;
      
      % Add father relationship to edges
    
      Mesh.Father_Edge(Mesh.Vert2Edge(i1,j1),1) = j1-nCoordinates;
      Mesh.Father_Edge(Mesh.Vert2Edge(j1,i2),1) = j1-nCoordinates;
      Mesh.Father_Edge(Mesh.Vert2Edge(i2,j2),1) = j2-nCoordinates;
      Mesh.Father_Edge(Mesh.Vert2Edge(j2,i3),1) = j2-nCoordinates;
      Mesh.Father_Edge(Mesh.Vert2Edge(i3,j3),1) = j3-nCoordinates;
      Mesh.Father_Edge(Mesh.Vert2Edge(j3,i4),1) = j3-nCoordinates;
      Mesh.Father_Edge(Mesh.Vert2Edge(i4,j4),1) = j4-nCoordinates;
      Mesh.Father_Edge(Mesh.Vert2Edge(j4,i1),1) = j4-nCoordinates;
      Mesh.Father_Edge(Mesh.Vert2Edge(j1,jc),2) = i;
      Mesh.Father_Edge(Mesh.Vert2Edge(j2,jc),2) = i;
      Mesh.Father_Edge(Mesh.Vert2Edge(j3,jc),2) = i;
      Mesh.Father_Edge(Mesh.Vert2Edge(j4,jc),2) = i;
    end
  end
    
return