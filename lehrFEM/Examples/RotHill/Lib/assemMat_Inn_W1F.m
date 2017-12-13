function varargout = assemMat_Inn_W1F(Mesh,EHandle,varargin)
%   
%   assembling routine for interior facet integrals and W1F-elements
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    EDGES        N-by-2 matrix specifying all edges of the mesh.
%    BDFLAGS      P-by-1 matrix specifying the boundary type of each
%                 boundary edge in the mesh.
%    VERT2EDGE    M-by-M sparse matrix which specifies whether the two
%                 vertices i and j are connected by an edge with number
%                 VERT2EDGE(i,j).
%    EDGE2ELEM    P-by-2 matrix connecting edges to elements. The first
%                 column specifies the left hand side element where the
%                 second column specifies the right hand side element.
%    EDGELOC      P-by-3 matrix connecting egdes to local edges of elements. 
%    NORMALS      P-by-2 matrix specifying the normals on each edge. The
%                 normals on interior edges are chosen such that they point
%                 from the element with the lower number to the element
%                 with the higher number and on boundary edges such that
%                 they point outside the domain.
%    MATCH        P-by-2 matrix specifying wheter the edge orientation of
%                 the current edge matches the orientation of the left and
%                 right hand side element. 
%

%   Copyright 2006-2009 Patrick Meury, Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nEdges = size(Mesh.Edges,1);        % Number of edges 
  nElements = size(Mesh.Elements,1);  % Number of elements
  
  % Allocate memory
  
  I = zeros(36*nEdges,1);
  J = zeros(36*nEdges,1);
  A = zeros(36*nEdges,1);

  % Check for element and boundary flags
  
  if(isfield(Mesh,'ElemFlag')),
    ElemFlag = Mesh.ElemFlag; 
  else
    ElemFlag = zeros(nElements,1);
  end
  BdFlags = Mesh.BdFlags;
  
  % Assemble element contributions
  
  loc = 1:36;
  last = 0;
  for i = 1:nEdges
     
    % Check for interior edge  
      
    if(BdFlags(i) >= 0)
       
      % Extract edge data
      
      Edge = Mesh.Coordinates(Mesh.Edges(i,:),:);
      Normal = Mesh.Normals(i,:);
      
      % Extract left and right hand side element data
      
      LData.Element = Mesh.Edge2Elem(i,1);
      LData.ElemFlag = ElemFlag(LData.Element);
      LData.Vertices = Mesh.Coordinates(Mesh.Elements(LData.Element,:),:);
      LData.EdgeLoc = Mesh.EdgeLoc(i,1);
      LData.Match = Mesh.EdgeOrient(i,1);
      
      RData.Element = Mesh.Edge2Elem(i,2);
      RData.ElemFlag = ElemFlag(RData.Element);
      RData.Vertices = Mesh.Coordinates(Mesh.Elements(RData.Element,:),:);
      RData.EdgeLoc = Mesh.EdgeLoc(i,2);
      RData.Match = Mesh.EdgeOrient(i,2);
      
      % Compute element contributions
      
      Aloc = EHandle(Edge,Normal,BdFlags(i), ...
                     LData,RData,varargin{:});
      
      % Add contributions to stiffness matrix
      idx_l = [Mesh.Vert2Edge(Mesh.Elements(LData.Element,2),Mesh.Elements(LData.Element,3)) ...
            Mesh.Vert2Edge(Mesh.Elements(LData.Element,3),Mesh.Elements(LData.Element,1)) ...
            Mesh.Vert2Edge(Mesh.Elements(LData.Element,1),Mesh.Elements(LData.Element,2))];
      
      idx_r = [Mesh.Vert2Edge(Mesh.Elements(RData.Element,2),Mesh.Elements(RData.Element,3)) ...
            Mesh.Vert2Edge(Mesh.Elements(RData.Element,3),Mesh.Elements(RData.Element,1)) ...
            Mesh.Vert2Edge(Mesh.Elements(RData.Element,1),Mesh.Elements(RData.Element,2))];
      
      idx = [idx_l idx_r];
      
      % Determine the orientation

      if(Mesh.Edges(idx_l(1),1)==Mesh.Elements(LData.Element,2)),  p1_l = 1;  else    p1_l = -1;  end
      if(Mesh.Edges(idx_l(2),1)==Mesh.Elements(LData.Element,3)),  p2_l = 1;  else    p2_l = -1;  end
      if(Mesh.Edges(idx_l(3),1)==Mesh.Elements(LData.Element,1)),  p3_l = 1;  else    p3_l = -1;  end
      
      if(Mesh.Edges(idx_r(1),1)==Mesh.Elements(RData.Element,2)),  p1_r = 1;  else    p1_r = -1;  end
      if(Mesh.Edges(idx_r(2),1)==Mesh.Elements(RData.Element,3)),  p2_r = 1;  else    p2_r = -1;  end
      if(Mesh.Edges(idx_r(3),1)==Mesh.Elements(RData.Element,1)),  p3_r = 1;  else    p3_r = -1;  end

      Peori = diag([p1_l p2_l p3_l p1_r p2_r p3_r]); % scaling matrix taking into account orientations
      Aloc = Peori*Aloc*Peori;
      Aloc = Aloc;
      
      I(loc) = set_Rows(idx,6);
      J(loc) = set_Cols(idx,6);
      A(loc) = Aloc(:);
      
      loc = loc+36;
      last = last+36;
      
    end
      
  end
  
  % Assign output arguments
  
  loc = 1:last;
  if(nargout > 1)
    varargout{1} = I(loc);
    varargout{2} = J(loc);
    varargout{3} = A(loc);
  else
    varargout{1} = sparse(I(loc),J(loc),A(loc));  
  end
  
return