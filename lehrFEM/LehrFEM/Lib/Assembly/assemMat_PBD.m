function varargout = assemMat_PBD(Mesh,EHandle,varargin)
% ASSEMMAT_LFE Assemble quadratic FE contributions.
%
%   A = ASSEMMAT_PBD(MESH,EHANDLE) assembles the global matrix from the
%   local element contributions given by the function handle EHANDLE and
%   returns the matrix in a sparse representation.
%
%   A = ASSEMMAT_PBD(MESH,EHANDLE,EPARAM) handles the variable length
%   argument list EPARAM to the function handle EHANDLE during the assembly
%   process. 
%
%   [I,J,A] = ASSEMMAT_PBD(MESH,EHANDLE) assembles the global matrix from
%   the local element contributions given by the function handle EHANDLE
%   and returns the matrix in an array representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%    EDGES        P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE    M-by-M sparse matrix which specifies wheter the two
%                 vertices i and j are connected by an edge with number
%                 VERT2EDGE(i,j).
%    EDGE2ELEM    N-by-2 matrix connecting edges to elements. The first
%                 column specifies the left hand side element where the
%                 second column specifies the right hand side element.
%    EDGELOC      P-by-2 matrix connecting egdes to local edges of
%                 elements.
%    DELTA        P-by-1 matrix specifying the boundary correction term on
%                 every edge.
%
%   Example:
%
%   A = assemMat_PBD(Mesh,@STIMA_Lapl_PBD,P7O6());
%  
%   See also set_Rows, set_Cols.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  
  % Preallocate memory
  
  I = zeros(36*nElements,1);
  J = zeros(36*nElements,1);
  A = zeros(36*nElements,1);
  
  % Assemble element contributions
  
  loc = 1:36;
  for i = 1:nElements
    
    % Extract vertices of current element
    
    vidx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(vidx,:);
      
    % Compute element contributions
    
    eidx = Mesh.Vert2Edge(vidx(2),vidx(3));
    Aloc = EHandle(Vertices,Mesh.ElemFlag(i),Mesh.Delta(eidx),varargin{:});
   
    % Add contributions to stiffness matrix
    
    idx = [vidx ...
           Mesh.Vert2Edge(vidx(1),vidx(2))+nCoordinates ...
           Mesh.Vert2Edge(vidx(2),vidx(3))+nCoordinates ...
           Mesh.Vert2Edge(vidx(3),vidx(1))+nCoordinates];
    I(loc) = set_Rows(idx,6);
    J(loc) = set_Cols(idx,6);
    A(loc) = Aloc(:);
    loc = loc+36;
    
  end
  
  % Assign output arguments
  
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
  else
    varargout{1} = sparse(I,J,A);      
  end
  
return