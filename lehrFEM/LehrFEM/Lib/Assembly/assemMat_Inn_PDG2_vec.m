function varargout = assemMat_Inn_PDG2_vec(Mesh,dim,EHandle,varargin)
%ASSEMMAT_INN_PDG2_VEC Assemble vector-valued DG interior edge contributions
%
%   A = ASSEMMAT_INN_PDG_VEC(MESH,DIM,EHANDLE) assembles the global matrix
%   from the local edge contributions fiven by the function handle EHANDLE
%   and returns the matrix in sparse representation.  This code is
%   for vector-valued functions, for example mixed DG, with the same basis
%   functions in each component of the image.  The structure of the
%   discrete vectors should be such that all coefficients corresponding to
%   the first dimension come first, then the second dimension, and so on.
%
%   A = ASSEMMAT_INN_PDG2_VEC(MESH,DIM,EHANDLE,EPARAM) passes the variable-
%   length argument list EPARAM to the function handle EHANDLE during the assembly
%   process.
%
%   [I,J,A] = ASSEMMAT_INN_PDG2_VEC(...) returns the matrix in an array
%   representation.
%
%   DIM is the dimension of the image space.  For example, for mixed DG,
%   the image space consists of one dimension for the function value and
%   two for the gradient, amounting to a total of three.  If the matrix is
%   not square, use DIM = [DIM1,DIM2], where DIM1 is the dimension of the
%   elements of the space determining the number of rows of the matrix (ie.
%   the dimension of the test functions in finite elements) and DIM2 is the
%   dimension of the elements of the space determining the number of
%   columns of the matrix (ie. the dimension of the solution).
%
%   EHANDLE must take at least the arguments:
%     EDGE        Coordinates of edge
%     NORMAL      Normal vector
%     EDGEDATA    Struct containing arbitrary edge information
%     LDATA       Struct containing information on left element
%     RDATA       Struct containing information on right element
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 or N-by-4 matrix specifying the elements of the 
%                 mesh.
%    EDGES        L-by-2 matrix specifying all edges of the mesh.
%    BDFLAGS      L-by-1 matrix specifying the boundary type of each
%                 boundary edge in the mesh.
%    VERT2EDGE    M-by-M sparse matrix which specifies whether the two
%                 vertices i and j are connected by an edge with number
%                 VERT2EDGE(i,j).
%    EDGE2ELEM    L-by-2 matrix connecting edges to elements. The first
%                 column specifies the left hand side element where the
%                 second column specifies the right hand side element.
%    EDGELOC      L-by-3 or P-by-4 matrix connecting egdes to local edges
%                 of elements. 
%    NORMALS      L-by-2 matrix specifying the normals on each edge. The
%                 normals on interior edges are chosen such that they point
%                 from the element with the lower number to the element
%                 with the higher number and on boundary edges such that
%                 they point outside the domain.
%    MATCH        L-by-2 matrix specifying whether the edge orientation of
%                 the current edge matches the orientation of the left and
%                 right hand side element.
%    ELEMDATA     N-by-1 structure array containing at least the fields:
%       NDOFS       The number of degrees of freedom on the corresponding
%                   element.
%    EDGEDATA     L-by-1 structure array containing at least the fields:
%       NDOFS       The total number of degrees of freedom on the two
%                   elements adjacent to a given edge.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  nEdges = size(Mesh.Edges,1);      % Number of edges
  dim1 = dim(1);
  dim2 = dim(end);
  
  % Find relevant edges
  doEdge = Mesh.BdFlags>=0;

  % Count degrees of freedom on all elements and calculate the number of 
  % matrix entries.
  nDofs = [Mesh.EdgeData.nDofs];     % Number of degrees of freedom on edge        
  nDofsSum = cumsum([0,[Mesh.ElemData.nDofs]]); % Number of degrees of freedom on all preceding elements
  nDofsTot = nDofsSum(end);          % Total number of degrees of freedom per dimension
  numel = dim1*dim2*sum(nDofs(doEdge).^2); % Number of matrix entries
  
  % Preallocate memory
  I = zeros(numel,1);
  J = zeros(numel,1);
  A = zeros(numel,1);
  
  % Assemble element contributions
  loc = 1;
  dim1Offset = nDofsTot*(0:dim1-1);
  dim2Offset = nDofsTot*(0:dim2-1);
  for i=1:nEdges
      
    % Check for interior edge  
    if(doEdge(i))
       
      % Extract edge data
      Edge = Mesh.Coordinates(Mesh.Edges(i,:),:);
      Normal = Mesh.Normals(i,:);
      
      % Extract left and right hand side element data
      LData.Element = Mesh.Edge2Elem(i,1);
      LData.ElemData = Mesh.ElemData(LData.Element);
      LData.Vertices = Mesh.Coordinates(Mesh.Elements(LData.Element,:),:);
      LData.EdgeLoc = Mesh.EdgeLoc(i,1);
      LData.Match = Mesh.EdgeOrient(i,1);
      
      RData.Element = Mesh.Edge2Elem(i,2);
      RData.ElemData = Mesh.ElemData(RData.Element);
      RData.Vertices = Mesh.Coordinates(Mesh.Elements(RData.Element,:),:);
      RData.EdgeLoc = Mesh.EdgeLoc(i,2);
      RData.Match = Mesh.EdgeOrient(i,2);
      
      % Compute element contributions
      Aloc = EHandle(Edge,Normal,Mesh.EdgeData(i), ...
                     LData,RData,varargin{:});
      
      % Add contributions to stiffness matrix
      idx_l = nDofsSum(LData.Element) + (1:LData.ElemData.nDofs); 
      idx_r = nDofsSum(RData.Element) + (1:RData.ElemData.nDofs);
      idx0 = [idx_l idx_r];
      idx1 = zeros(1,dim1*nDofs(i));
      idx1(:) = idx0(ones(dim1,1),:)' + dim1Offset(ones(nDofs(i),1),:);
      idx2 = zeros(1,dim2*nDofs(i));
      idx2(:) = idx0(ones(dim2,1),:)' + dim2Offset(ones(nDofs(i),1),:);
      loc1 = loc + dim1*dim2*nDofs(i)^2 - 1;
  
      I(loc:loc1) = set_Rows(idx1,dim2*nDofs(i));
      J(loc:loc1) = set_Cols(idx2,dim1*nDofs(i));
      A(loc:loc1) = Aloc(:);
      
      loc = loc1 + 1;
      
    end
      
  end
  
  % Assign output arguments
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
  else
    varargout{1} = sparse(I,J,A,nDofsSum(end),nDofsSum(end));  
  end
  
return