function varargout = assemMat_Inn_PDG(Mesh,nDofs,EHandle,varargin)
% ASSEMMAT_INN_PDG Assemble interior edge contributions.
%
%   A = ASSEMMAT_INN_PDG(MESH,NDOFS,EHANDLE) assembles the global matrix
%   from the local element contributions given by the function handle
%   EHANDLE and returns the matrix in a sparse representation. The integer
%   NDOFS specifies the number of degrees of freedom per element.
%
%   A = ASSEMMAT_INN_PDG(MESH,NDOFS,EHANDLE,EPARAM) handles the variable
%   length argument list EPARAM to the function handle EHANDLE during the
%   assembly process. 
%
%   [I,J,A] = ASSEMMAT_INN_PDG(MESH,NDOFS,EHANDLE) assembles the global
%   matrix from the local element contributions given by the function
%   handle EHANDLE and returns the matrix in an array representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 or N-by-4 matrix specifying the elements of the 
%                 mesh.
%    EDGES        N-by-2 matrix specifying all edges of the mesh.
%    BDFLAGS      P-by-1 matrix specifying the boundary type of each
%                 boundary edge in the mesh.
%    VERT2EDGE    M-by-M sparse matrix which specifies whether the two
%                 vertices i and j are connected by an edge with number
%                 VERT2EDGE(i,j).
%    EDGE2ELEM    P-by-2 matrix connecting edges to elements. The first
%                 column specifies the left hand side element where the
%                 second column specifies the right hand side element.
%    EDGELOC      P-by-3 or P-by-4 matrix connecting egdes to local edges
%                 of elements. 
%    NORMALS      P-by-2 matrix specifying the normals on each edge. The
%                 normals on interior edges are chosen such that they point
%                 from the element with the lower number to the element
%                 with the higher number and on boundary edges such that
%                 they point outside the domain.
%    MATCH        P-by-2 matrix specifying whether the edge orientation of
%                 the current edge matches the orientation of the left and
%                 right hand side element. 
%
%   Example:
%
%   A = assemMat_Inn_PDG(Mesh,3,@STIMA_Lapl_Inn_PDG,1,P3O3(), ...
%                        @shap_DGCR,@grad_shap_DGCR);
%  
%   See also set_Rows, set_Cols.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nEdges = size(Mesh.Edges,1);        % Number of edges 
  nElements = size(Mesh.Elements,1);  % Number of elements
  
  % Allocate memory
  
  I = zeros((2*nDofs)^2*nEdges,1);
  J = zeros((2*nDofs)^2*nEdges,1);
  A = zeros((2*nDofs)^2*nEdges,1);

  % Check for element and boundary flags
  
  if(isfield(Mesh,'ElemFlag')),
    ElemFlag = Mesh.ElemFlag; 
  else
    ElemFlag = zeros(nElements,1);
  end
  BdFlags = Mesh.BdFlags;
  
  % Assemble element contributions
  
  loc = 1:((2*nDofs)^2);
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
      
      
      
      idx_l = nDofs*(LData.Element-1) + (1:nDofs); 
      idx_r = nDofs*(RData.Element-1) + (1:nDofs);
      idx = [idx_l idx_r];
  
            
      I(loc) = set_Rows(idx,2*nDofs);
      J(loc) = set_Cols(idx,2*nDofs);
      A(loc) = Aloc(:);
      
      loc = loc + (2*nDofs)^2;
      last = last + (2*nDofs)^2;
      
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