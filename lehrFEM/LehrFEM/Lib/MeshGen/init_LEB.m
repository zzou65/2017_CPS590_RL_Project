function New_Mesh = init_LEB(Old_Mesh)
% INIT_LEB Initialize largest edge bisection.
%  
%   MESH = INIT_LEB(MESH) initializes the struct MESH for largest edge
%   bisections.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh. 
%    EDGES       P-by-2 matrix specifying the edges of the mesh. 
%    BDFLAGS     P-by-1 matrix specifying additional edge information for each
%                edge of the mesh.
%    VERT2EDGE   M-by-M sparse matrix which specifies wheter the two vertices
%                i and j are connected by an edge with number VERT2EDGE(i,j).
%
%   Example:
%
%   Mesh = init_LEB(Mesh);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  nElements = size(Old_Mesh.Elements,1);
  nEdges = size(Old_Mesh.Edges,1);
 
  % Copy data from old mesh to new mesh

  New_Mesh.Coordinates = Old_Mesh.Coordinates;
  New_Mesh.Elements = Old_Mesh.Elements;
  if(isfield(Old_Mesh,'ElemFlag'))
    New_Mesh.ElemFlag = Old_Mesh.ElemFlag;  
  end
  if(isfield(Old_Mesh,'Edges') && isfield(Old_Mesh,'Vert2Edge'))
    New_Mesh.Edges = Old_Mesh.Edges;
    New_Mesh.Vert2Edge = Old_Mesh.Vert2Edge;
  end
  if(isfield(Old_Mesh,'BdFlags'))
    New_Mesh.BdFlags = Old_Mesh.BdFlags;  
  end
  
  % Preallocate memory
  
  New_Mesh.Neigh = zeros(nElements,3);
  New_Mesh.Opp_Vert = zeros(nElements,3);
  Edge2Elem = zeros(nEdges,2);
  LR_Vert = zeros(nEdges,2);
    
  for i = 1:nElements     
    id = New_Mesh.Elements(i,:);
    
    % Reorder current element such that edge 3 is the longest edge
    
    [h_max,id_max] = max([ norm(Old_Mesh.Coordinates(id(3),:)-Old_Mesh.Coordinates(id(2),:)), ...
                           norm(Old_Mesh.Coordinates(id(1),:)-Old_Mesh.Coordinates(id(3),:)), ...
                           norm(Old_Mesh.Coordinates(id(2),:)-Old_Mesh.Coordinates(id(1),:)) ]);
    id = id([rem(id_max+3,3)+1 rem(id_max+4,3)+1 id_max]);
    New_Mesh.Elements(i,:) = id;
    
    % Build auxilliary connectivity tables
    
    j1 = Old_Mesh.Vert2Edge(id(2),id(3));
    if(id(2) < id(3))
      Edge2Elem(j1,1) = i;
      Opp_Vert(j1,1) = 1;
    else
      Edge2Elem(j1,2) = i;
      Opp_Vert(j1,2) = 1;
    end
    j2 = Old_Mesh.Vert2Edge(id(3),id(1));
    if(id(3) < id(1))
      Edge2Elem(j2,1) = i; 
      Opp_Vert(j2,1) = 2;
    else
      Edge2Elem(j2,2) = i; 
      Opp_Vert(j2,2) = 2;
    end
    j3 = Old_Mesh.Vert2Edge(id(1),id(2));
    if(id(1) < id(2))
      Edge2Elem(j3,1) = i;
      Opp_Vert(j3,1) = 3;
    else
      Edge2Elem(j3,2) = i;
      Opp_Vert(j3,2) = 3;
    end
  end
  
  % Connect elements to neighbours and opposite vertices
  
  for i = 1:nEdges
    Elem_Left = Edge2Elem(i,1);
    Elem_Right = Edge2Elem(i,2);
    if(Elem_Left > 0 && Elem_Right > 0)
      New_Mesh.Neigh(Elem_Left,Opp_Vert(i,1)) = Elem_Right;
      New_Mesh.Opp_Vert(Elem_Left,Opp_Vert(i,1)) = Opp_Vert(i,2);
      New_Mesh.Neigh(Elem_Right,Opp_Vert(i,2)) = Elem_Left;
      New_Mesh.Opp_Vert(Elem_Right,Opp_Vert(i,2)) = Opp_Vert(i,1);
    elseif(Elem_Left > 0)
      New_Mesh.Neigh(Elem_Left,Opp_Vert(i,1)) = Old_Mesh.BdFlags(i);
    else
      New_Mesh.Neigh(Elem_Right,Opp_Vert(i,2)) = Old_Mesh.BdFlags(i); 
    end
  end
  
return
