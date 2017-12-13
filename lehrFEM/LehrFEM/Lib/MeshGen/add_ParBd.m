function Mesh = add_ParBd(Mesh,DHandle,varargin)
% ADD_PARBD Add parabolic boundaries to the mesh.
%
%   MESH = ADD_PARBD(MESH,DHANDLE) adds parabolic boundaries to the struct
%   MESH.
%
%   MESH = ADD_PARBD(MESH,DHANDLE,DPARAM) also handles the variable length
%   argument list DPARAM to the signed distance function DHANDLE.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh.   
%    EDGES       P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE   M-by-M sparse matrix which specifies wheter the two
%                vertices i and j are connected by an edge with number
%                VERT2EDGE(i,j).
%    EDGE2ELEM   N-by-2 matrix connecting edges to elements. The first
%                column specifies the left hand side element where the
%                second column specifies the right hand side element.
%    EDGELOC     P-by-2 matrix connecting egdes to local edges of elements. 
%   
%   Parabloic boundaries can only be appended to mesh where all elements
%   share at most one edge with the boundary.
%
%   Example:
%
%   Mesh = add_ParBd(Mesh,@dist_circ,[0 0],1);
%
%   See also get_BdEdges, add_Edge2Elem.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);

  % Check number of boundary edge 
  
  BdEdges = zeros(nElements,1);
  Loc = get_BdEdges(Mesh);
  ElemLeft = Mesh.Edge2Elem(Loc,1);
  ElemLeft = ElemLeft(ElemLeft > 0);
  BdEdges(ElemLeft) = 1;
  ElemRight = Mesh.Edge2Elem(Loc,2);
  ElemRight = ElemRight(ElemRight > 0);
  BdEdges(ElemRight) = 1;
  if(max(BdEdges) > 1)
    error('More then one boundary edge per element');  
  end
  
  % Compute distance parameter
  
  Mesh.Delta = zeros(nEdges,1);
  for i = 1:nElements      
    Elem = Mesh.Elements(i,:);    
    if(BdEdges(i) > 0)
      while(Mesh.BdFlags(Mesh.Vert2Edge(Elem(2),Elem(3))) >= 0)
        tmp = Elem(1);
        Elem(1) = Elem(2);
        Elem(2) = Elem(3);
        Elem(3) = tmp;
      end
      Mesh.Elements(i,:) = Elem;
      x = (Mesh.Coordinates(Elem(2),:)+Mesh.Coordinates(Elem(3),:))/2;
      Mesh.Delta(Mesh.Vert2Edge(Elem(2),Elem(3))) = -DHandle(x,varargin{:});
    else
      Mesh.Delta(Mesh.Vert2Edge(Elem(2),Elem(3))) = 0;
    end
  end
  
  % Update mesh data structure
  
  Mesh = add_Edge2Elem(Mesh);
  
return