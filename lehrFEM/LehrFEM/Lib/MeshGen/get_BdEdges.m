function Loc = get_BdEdges(Mesh)
% GET_BDEDGES Extract boundary edges of the mesh.
%
%   LOC = GET_BDEDGES(MESH) Extract boundary edge locations from the mesh.
%
%   The struct mesh should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the mesh. 
%    EDGES       P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE   M-by-M sparse matrix which specifies whether the two vertices
%                i and j are connected by an edge with number VERT2EDGE(i,j).
%
%   Example:
%
%   loc = get_BdEdges(Mesh);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  nCoordinates = size(Mesh.Coordinates,1);
  [nElements,nVert] = size(Mesh.Elements);
  nEdges = size(Mesh.Edges,1);

  if(isfield(Mesh,'BdFlags'))
    
    % Boundary information is available  
      
    Loc = find(Mesh.BdFlags < 0);
  else
    
    % No boundary information available  
      
    EdgeMarkers = zeros(nEdges,1);
    if(nVert == 3)
        
      % Triangular elements  
        
      for i = 1:nElements
        i1 = Mesh.Elements(i,1);
        i2 = Mesh.Elements(i,2);
        i3 = Mesh.Elements(i,3);
        EdgeMarkers(Mesh.Vert2Edge(i1,i2)) = EdgeMarkers(Mesh.Vert2Edge(i1,i2))+1;
        EdgeMarkers(Mesh.Vert2Edge(i2,i3)) = EdgeMarkers(Mesh.Vert2Edge(i2,i3))+1;
        EdgeMarkers(Mesh.Vert2Edge(i3,i1)) = EdgeMarkers(Mesh.Vert2Edge(i3,i1))+1;
      end
      Loc = find(EdgeMarkers == 1);
    else
       
      % Quadrilateral elements  
        
      for i = 1:nElements
        i1 = Mesh.Elements(i,1);
        i2 = Mesh.Elements(i,2);
        i3 = Mesh.Elements(i,3);
        i4 = Mesh.Elements(i,4);
        EdgeMarkers(Mesh.Vert2Edge(i1,i2)) = EdgeMarkers(Mesh.Vert2Edge(i1,i2))+1;
        EdgeMarkers(Mesh.Vert2Edge(i2,i3)) = EdgeMarkers(Mesh.Vert2Edge(i2,i3))+1;
        EdgeMarkers(Mesh.Vert2Edge(i3,i4)) = EdgeMarkers(Mesh.Vert2Edge(i3,i4))+1;
        EdgeMarkers(Mesh.Vert2Edge(i4,i1)) = EdgeMarkers(Mesh.Vert2Edge(i4,i1))+1;
      end
    end
     Loc = find(EdgeMarkers == 1);
  end
       
return