function New_Mesh = add_VBdFlags(Old_Mesh,inflow_flag, outflow_flag)
% add_VBdFlags add array containing the boundary flags for vertices
%
%   MESH = add_VBdFlags(MESH) add array containing the boundary flags for vertices 
%  The struct MESH must at least contain the following fields:
%  
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh. 
%    EDGES       N-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE   M-by-M sparse matrix which specifies whether the two
%                vertices i and j are connected by an edge with number
%                VERT2EDGE(i,j).
%    EDGE2ELEM   P-by-2 matrix connecting edges to elements. The first
%                column specifies the left hand side element where the
%                second column specifies the right hand side element.
%   
%   Example:
%
%   Mesh = add_VBdFlags(Mesh,v);

%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  nCoordinates = size(Old_Mesh.Coordinates,1);
  nElements=size(Old_Mesh.Elements,1);
  % Copy old mesh
  
  New_Mesh = Old_Mesh;
  
  % Allocate memory for edge normals
  
  New_Mesh.VBdFlags = zeros(nCoordinates,1);
  
  BndEdges=get_BdEdges(Old_Mesh);
  for i = BndEdges'
      
    vid=New_Mesh.Edges(i,:); 
    New_Mesh.VBdFlags(vid(1)) = New_Mesh.VBdFlags(vid(1)) + New_Mesh.BdFlags(i);
    New_Mesh.VBdFlags(vid(2)) = New_Mesh.VBdFlags(vid(2)) + New_Mesh.BdFlags(i);    
  end
  
return