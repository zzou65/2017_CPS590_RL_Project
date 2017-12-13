function New_Mesh = add_Edge2Elem(Old_Mesh)
% ADD_EDGE2ELEM Connects edges to elements.
%
%   MESH = ADD_EDGE2ELEM(MESH) adds an additional connectivity table to the
%   struct MESH:
%    EDGE2ELEM P-by-2 matrix connecting edges to elements. The first column
%              specifies the left hand side element and the second column
%              specifies the right hand side element.
%    EDGELOC   P-by-3 or P-by-4 matrix connecting egdes to local edges of
%              elements. 
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the mesh.   
%    EDGES       P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE   M-by-M sparse matrix which specifies wheter the two vertices i
%                and j are connected by an edge with number VERT2EDGE(i,j).
%
%   Note that the elements of MESH should be oriented counterclockwise.
%
%   Example:
%
%   Mesh = add_Edge2Elem(Mesh);
%
%   See also add_Edges, orient_Elems.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  [nElements,nVert] = size(Old_Mesh.Elements);
  nEdges = size(Old_Mesh.Edges,1);
  
  % Preallocate memory
  
  Edge2Elem = zeros(nEdges,2);
  EdgeLoc = zeros(nEdges,2);
  
  % Connect edges to elements
  
  if(nVert == 3)
      
    % Triangular elements
      
    for i = 1:nElements
      idx_1 = Old_Mesh.Elements(i,1);
      idx_2 = Old_Mesh.Elements(i,2);
      idx_3 = Old_Mesh.Elements(i,3);
      EdgeId_1 = Old_Mesh.Vert2Edge(idx_2,idx_3);
      EdgeId_2 = Old_Mesh.Vert2Edge(idx_3,idx_1);
      EdgeId_3 = Old_Mesh.Vert2Edge(idx_1,idx_2);
      %[i, idx_1, idx_2, idx_3]
      if(idx_2 < idx_3)
      %  disp([num2str(i) ' is left of (' num2str(idx_2) ',' num2str(idx_3) ')']);
        Edge2Elem(EdgeId_1,1) = i;
        EdgeLoc(EdgeId_1,1) = 1;
      else
      %  disp([num2str(i) ' is right of (' num2str(idx_3) ',' num2str(idx_2) ')']);
        Edge2Elem(EdgeId_1,2) = i;
        EdgeLoc(EdgeId_1,2) = 1;
      end
      if(idx_3 < idx_1)
      %  disp([num2str(i) ' is left of (' num2str(idx_3) ',' num2str(idx_1) ')']);
        Edge2Elem(EdgeId_2,1) = i;
        EdgeLoc(EdgeId_2,1) = 2;
      else
      %  disp([num2str(i) ' is right of (' num2str(idx_1) ',' num2str(idx_3) ')']);
        Edge2Elem(EdgeId_2,2) = i;
        EdgeLoc(EdgeId_2,2) = 2;
      end
      if(idx_1 < idx_2)
      %  disp([num2str(i) ' is left of (' num2str(idx_1) ',' num2str(idx_2) ')']);
        Edge2Elem(EdgeId_3,1) = i;
        EdgeLoc(EdgeId_3,1) = 3;
      else
      %  disp([num2str(i) ' is right of (' num2str(idx_2) ',' num2str(idx_1) ')']);
        Edge2Elem(EdgeId_3,2) = i;
        EdgeLoc(EdgeId_3,2) = 3;
      end     
    end
      
    % added, making sure that bdry edge always occure in the first entry
    for k = 1:nEdges
      if(Edge2Elem(k,1)==0)
          Edge2Elem(k,1) = Edge2Elem(k,2);
          EdgeLoc(k,1) = EdgeLoc(k,2);
          Edge2Elem(k,2) = 0;
          EdgeLoc(k,2) = 0;
      end
    end
    
  else
      
    % Quadrilateral elements 
    for i = 1:nElements
      idx_1 = Old_Mesh.Elements(i,1);
      idx_2 = Old_Mesh.Elements(i,2);
      idx_3 = Old_Mesh.Elements(i,3);
      idx_4 = Old_Mesh.Elements(i,4);
      EdgeId_1 = Old_Mesh.Vert2Edge(idx_1,idx_2);
      EdgeId_2 = Old_Mesh.Vert2Edge(idx_2,idx_3);
      EdgeId_3 = Old_Mesh.Vert2Edge(idx_3,idx_4);
      EdgeId_4 = Old_Mesh.Vert2Edge(idx_4,idx_1);
      if(idx_1 < idx_2)
        Edge2Elem(EdgeId_1,1) = i;
        EdgeLoc(EdgeId_1,1) = 1;
      else
        Edge2Elem(EdgeId_1,2) = i;
        EdgeLoc(EdgeId_1,2) = 1;
      end
      if(idx_2 < idx_3)
        Edge2Elem(EdgeId_2,1) = i;
        EdgeLoc(EdgeId_2,1) = 2;
      else
        Edge2Elem(EdgeId_2,2) = i;
        EdgeLoc(EdgeId_2,2) = 2;
      end
      if(idx_3 < idx_4)
        Edge2Elem(EdgeId_3,1) = i;
        EdgeLoc(EdgeId_3,1) = 3;
      else
        Edge2Elem(EdgeId_3,2) = i;
        EdgeLoc(EdgeId_3,2) = 3;
      end
      if(idx_4 < idx_1)
        Edge2Elem(EdgeId_4,1) = i;
        EdgeLoc(EdgeId_4,1) = 4;   
      else
        Edge2Elem(EdgeId_4,2) = i;
        EdgeLoc(EdgeId_4,2) = 4;  
      end
    end
      
    % added, making sure that bdry edge always occure in the first entry
    for k = 1:nEdges
      if(Edge2Elem(k,1)==0)
          Edge2Elem(k,1) = Edge2Elem(k,2);
          EdgeLoc(k,1) = EdgeLoc(k,2);
          Edge2Elem(k,2) = 0;
          EdgeLoc(k,2) = 0;
      end
    end
    
  end
  
  % Assign output arguments
  
  New_Mesh = deal(Old_Mesh);
  New_Mesh.Edge2Elem = Edge2Elem;
  New_Mesh.EdgeLoc = EdgeLoc;
  
return
