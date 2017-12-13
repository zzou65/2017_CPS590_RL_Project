function New_Mesh = add_Edges(Old_Mesh)
% ADD_EDGES Add edge lists.
%
%   MESH = ADD_EDGES(MESH) adds additional edge information to the struct MESH:
%    EDGES     N-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE M-by-M sparse matrix which specifies whether the two vertices i
%              and j are connected by an edge with number VERT2EDGE(i,j).
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the mesh. 
%
%   Example:
%
%   Mesh = add_Edges(Mesh);

%   Copyright 2005-2007 Patrick Meury & Holger Heumann & Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  nCoordinates = size(Old_Mesh.Coordinates,1);
  [nElements,nVert] = size(Old_Mesh.Elements); 
  
  % Initialize constants
  
  if(isfield(Old_Mesh,'Max_Nodes'))
    MAX_NODES = Old_Mesh.Max_Nodes;
  else
    MAX_NODES = 6;
  end
  MAX_EDGES = nVert*nElements; 
 
  % Build list of adjacent nodes
  
  AdjNodes = zeros(nCoordinates,MAX_NODES);
  nAdjNodes = zeros(nCoordinates,1);
  if(nVert == 3)
      
    % Triangular elements  
      
    for i = 1:nElements
      for j = 1:3
        idx_1 = Old_Mesh.Elements(i,j);
        idx_2 = Old_Mesh.Elements(i,rem(j+3,3)+1);
        idx_3 = Old_Mesh.Elements(i,rem(j+4,3)+1);
      
        % Check for vertices 2 and 3 in the set of adjacent nodes of vertex 1
      
        sz = nAdjNodes(idx_1);
        tmp = AdjNodes(idx_1,:);
        occ_1 = 0;
        occ_2 = 0;
        k = 1;
        while((~occ_1 || ~occ_2) && k <= sz)
          occ_1 = max(occ_1,idx_2 == tmp(k));
          occ_2 = max(occ_2,idx_3 == tmp(k));
          k = k+1;
        end
        if(~occ_1)
          sz = sz+1;
          tmp(sz) = idx_2;    
        end
        if(~occ_2)
          sz = sz+1; 
          tmp(sz) = idx_3;
        end
        
        % Increase size of AdjNodes if tmp is too long
        
        if(sz > MAX_NODES)
          AdjNodes(idx_1,sz) = 0;
          MAX_NODES = sz;
        end
        
        % Update list of adjacent vertices
   
        nAdjNodes(idx_1) = sz;
        AdjNodes(idx_1,:) = tmp;
      end
    end
  else
      
    % Quadrilateral elements  
    
    for i = 1:nElements
      for j = 1:4
        idx_1 = Old_Mesh.Elements(i,j);
        idx_2 = Old_Mesh.Elements(i,rem(j+4,4)+1);
        idx_4 = Old_Mesh.Elements(i,rem(j+6,4)+1);
      
        % Check for vertices 2, 3 and 4 in the set of adjacent nodes of vertex 1
      
        sz = nAdjNodes(idx_1);
        tmp = AdjNodes(idx_1,:);
        occ_1 = 0;
        occ_2 = 0;
        k = 1;
        while((~occ_1 || ~occ_2) && k <= sz)
          occ_1 = max(occ_1,idx_2 == tmp(k));
          occ_2 = max(occ_2,idx_4 == tmp(k));
          k = k+1;
        end
        if(~occ_1)
          sz = sz+1;
          tmp(sz) = idx_2;    
        end
        if(~occ_2)
          sz = sz+1; 
          tmp(sz) = idx_4;
        end
        
        % Increase size of AdjNodes if tmp is too long
        
        if(sz > MAX_NODES)
          AdjNodes(idx_1,sz) = 0;
          MAX_NODES = sz;
        end
        
        % Update list of adjacent vertices
           
        nAdjNodes(idx_1) = sz;
        AdjNodes(idx_1,:) = tmp;
      end
    end
  end
  
  % Squeeze adjacency list
  
  MAX_NODES = max(nAdjNodes);
  
  AdjNodes = AdjNodes(:,1:MAX_NODES);
  
  % Build list of edges 
  
  nEdges = 0;
  I = zeros(1,MAX_EDGES);
  J = zeros(1,MAX_EDGES);
  EdgeNr = zeros(1,MAX_EDGES);
  Edges = zeros(MAX_EDGES,2);
  for i = 1:nCoordinates
    for j = 1:nAdjNodes(i)
      if(i < AdjNodes(i,j))
        nEdges = nEdges+1;
        Edges(nEdges,:) = [i AdjNodes(i,j)];
        I(nEdges) = i;
        J(nEdges) = AdjNodes(i,j);
        EdgeNr(nEdges) = nEdges;
      end
    end
  end
  
  % Squeeze arrays
  
  loc = 1:nEdges;
  I = I(loc);
  J = J(loc);
  EdgeNr = EdgeNr(loc);
  Edges = Edges(loc,:);
  
  % Add missing data and convert into sparse format
  
  EdgeNr = sparse([I J],[J I],[EdgeNr EdgeNr]);
      
  % Assign output arguments
  
  New_Mesh = deal(Old_Mesh);
  New_Mesh.Edges = Edges;
  New_Mesh.Vert2Edge = EdgeNr;
  New_Mesh.Max_Nodes=MAX_NODES;  
return
