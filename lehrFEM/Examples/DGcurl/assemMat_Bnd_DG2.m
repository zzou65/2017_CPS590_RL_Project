function varargout = assemMat_Bnd_DG2(Mesh,EHandle,varargin)
% ASSEMMAT_BND_DG2 Assemble vectorial discontinuous Lagrangian nodal FE boundary edge
% contributions.
%
%   A = ASSEMMAT_BND_DG2(MESH,EHANDLE) assembles the global matrix from the
%   local element contributions given by the function handle EHANDLE and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMAT_BND_DG(MESH,EHANDLE) assembles the global matrix
%   from the local element contributions given by the function handle
%   EHANDLE and returns the matrix in an array representation.
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

%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);  % Number of elements in the mesh
  nEdges = size(Mesh.Edges,1);        % Number of edges in the mesh

  % Allocate memory

  I = zeros(36*nEdges,1);
  J = zeros(36*nEdges,1);
  A = zeros(36*nEdges,1);
  
  % Check for element flags
  
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
  
    % Check for boundary edge  
      
    if(BdFlags(i) < 0)  
    
      Edge = Mesh.Coordinates(Mesh.Edges(i,:),:);
      Normal = Mesh.Normals(i,:);
        
      % Extract left or right hand side element data
    
      if(Mesh.Edge2Elem(i,1) > 0)
        Data.Element = Mesh.Edge2Elem(i,1);
        Data.ElemFlag = ElemFlag(Data.Element);
        Data.Vertices = Mesh.Coordinates(Mesh.Elements(Data.Element,:),:);
        Data.EdgeLoc = Mesh.EdgeLoc(i,1);
        Data.Match = Mesh.EdgeOrient(i,1);
      else
        Data.Element = Mesh.Edge2Elem(i,2);
        Data.ElemFlag = ElemFlag(Data.Element);
        Data.Vertices = Mesh.Coordinates(Mesh.Elements(Data.Element,:),:);
        Data.EdgeLoc = Mesh.EdgeLoc(i,2);
        Data.Match = Mesh.EdgeOrient(i,2);
      end
            
      % Compute element contributions
    
      Aloc = EHandle(Edge,Normal,BdFlags(i), ...
                     Data,varargin{:});   
    
      % Add element contributions to stiffness matrix
    
      idx = 6*(Data.Element-1)+[1 2 3 4 5 6];
    
      I(loc) = set_Rows(idx,6);
      J(loc) = set_Cols(idx,6);
      A(loc) = Aloc(:);
    
      loc = loc + 36;  
      last = last + 36;
      
    end
      
  end
  
  % Assign output arguments
  
  loc = 1:last;
  if(nargout > 1)
    varargout{1} = I(loc);
    varargout{2} = J(loc);
    varargout{3} = A(loc);
  else
    varargout{1} = sparse(I(loc),J(loc),A(loc),6*nElements,6*nElements);  
  end
  
return
