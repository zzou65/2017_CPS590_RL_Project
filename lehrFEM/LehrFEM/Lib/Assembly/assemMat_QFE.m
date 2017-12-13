function varargout = assemMat_QFE(Mesh,EHandle,varargin)
% ASSEMMAT_QFE Assemble quadratic FE contributions.
%
%   A = ASSEMMAT_QFE(MESH,EHANDLE) assembles the global matrix from the 
%   local element contributions given by the function handle EHANDLE and
%   returns the matrix in a sparse representation.
%
%   A = ASSEMMAT_QFE(MESH,EHANDLE,EPARAM) handles the variable length 
%   argument list EPARAM to the function handle EHANDLE during the assembly
%   process. 
%
%   [I,J,A] = ASSEMMAT_QFE(MESH,EHANDLE) assembles the global matrix
%   from the local element contributions given by the function handle
%   EHANDLE and returns the matrix in an array representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%    EDGES        P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE    M-by-M sparse matrix which specifies wheter the two
%                 vertices i and j are connected by an edge with number
%                 VERT2EDGE(i,j).
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   Mesh = add_Edges(Mesh);
%   Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
%   A = assemMat_QFE(Mesh,@STIMA_Lapl_QFE);
%
%   See also set_Rows, set_Cols.

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
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
     
    % Extract vertices of the current element
  
    vidx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(vidx,:);
    
    % Compute element contributions
    
    Aloc = EHandle(Vertices,Mesh.ElemFlag(i),varargin{:});
      
    % Extract global edge numbers
    
    idx = [vidx ...
           Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))+nCoordinates ...
           Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3))+nCoordinates ...
           Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1))+nCoordinates];
    
    % Add contributions to global matrix
    
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