function New_Mesh = add_DGData(Old_Mesh)
% ADD_NORMALS Add edge normals.
%
%   MESH = ADD_DGDATA(MESH) adds the normals and relative edge orientations
%   to the struct MESH:
%    NORMALS P-by-2 matrix specifying the normals on each edge. The normals
%            on interior edges are chosen such that they point from the
%            element with the higher number to the element with the lower
%            number and on boundary edges such that they point inside the
%            domain.
%            Note that, even if stated otherwise, this is the orientation
%            of the normal vectors that tends to work.
%    MATCH   P-by-2 matrix specifying whether the edge orientation of the
%            current edge matches the orientation of the left and right
%            hand side element.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the 
%                mesh. 
%    EDGES       N-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE   M-by-M sparse matrix which specifies whether the two
%                vertices i and j are connected by an edge with number
%                VERT2EDGE(i,j).
%    EDGE2ELEM   P-by-2 matrix connecting edges to elements. The first
%                column specifies the left hand side element where the
%                second column specifies the right hand side element.
%    EDGELOC     P-by-3 matrix connecting egdes to local edges of elements. 
%   
%   Example:
%
%   Mesh = add_DGData(Mesh);

%   Copyright 2006-2007 Patrick Meury & Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  nEdges = size(Old_Mesh.Edges,1);
  nVert = size(Old_Mesh.Elements,2);
  
  % Copy old mesh
  
  New_Mesh = Old_Mesh;
  
  % Allocate memory for edge normals
  
  New_Mesh.Normals = zeros(nEdges,2);
  New_Mesh.EdgeOrient = zeros(nEdges,2);
  
  % Generate edges according to edge orientation  
  
  Rot = [0 1; -1 0];
  for i = 1:nEdges
    
    % Extract start and end point of the current edge  
      
    id_s = New_Mesh.Edges(i,1);
    id_e = New_Mesh.Edges(i,2);
      
    % Compute left hand side unit normal and invert it since, due to a bug,
    % we actually need the right hand side unit normal.
    
    normal = New_Mesh.Coordinates(id_e,:) - ...
             New_Mesh.Coordinates(id_s,:);
    normal = -normal*Rot/norm(normal);
         
    % Check whether edge is located on the boundary
    
    ELeft = New_Mesh.Edge2Elem(i,1);
    ERight = New_Mesh.Edge2Elem(i,2);
    if(ELeft == 0 || ERight == 0)
      if(ERight == 0)
        normal = -normal;
      end
    else   
      if(ELeft < ERight)
        normal = -normal;  
      end
    end
    
    % Compute relative edge orientation wrt. left and right hand side
    % element
    
    if(nVert==3)  % triangular elements
  
      LMatch = 0;
      if(ELeft > 0)
        Lloc = New_Mesh.EdgeLoc(i,1);
        if(id_s == New_Mesh.Elements(ELeft,rem(Lloc,3)+1))
          LMatch = 1;   
        else
          LMatch = -1;
        end
      end

      RMatch = 0;
      if(ERight > 0)
        Rloc = New_Mesh.EdgeLoc(i,2);
        if(id_s == New_Mesh.Elements(ERight,rem(Rloc,3)+1))
          RMatch = 1;
        else
          RMatch = -1;  
        end
      end
      
    else % quadrilateral elements
      
      LMatch = 0;
      if(ELeft > 0)
        Lloc = New_Mesh.EdgeLoc(i,1);
        if(id_s == New_Mesh.Elements(ELeft,Lloc))
          LMatch = 1;   
        else
          LMatch = -1;
        end
      end

      RMatch = 0;
      if(ERight > 0)
        Rloc = New_Mesh.EdgeLoc(i,2);
        if(id_s == New_Mesh.Elements(ERight,Rloc))
          RMatch = 1;
        else
          RMatch = -1;  
        end
      end
      
    end
    
    % Assign outer unit normal and relative edge orientation
    
    New_Mesh.Normals(i,:) = normal;
    New_Mesh.EdgeOrient(i,1) = LMatch;
    New_Mesh.EdgeOrient(i,2) = RMatch;
        
  end
  
return