function varargout = assemMat_BFE(Mesh,EHandle,varargin)
% ASSEMMAT_BFE Assemble bi-linear FE contributions.
%
%   A = ASSEMMAT_BFE(MESH,EHANDLE) assembles the global matrix from the
%   local element contributions given by the function handle EHANDLE and
%   returns the matrix in a sparse representation.
%
%   A = ASSEMMAT_BFE(MESH,EHANDLE,EPARAM) handles the variable length
%   argument list EPARAM to the function handle EHANDLE during the assembly
%   process. 
%
%   [I,J,A] = ASSEMMAT_BFE(MESH,EHANDLE) assembles the global matrix from
%   the local element contributions given by the function handle EHANDLE
%   and returns the matrix in an array representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-4 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%   Example:
%
%   Mesh = load_Mesh('Coord_Sqr_QElem.dat','Elem_Sqr_QElem.dat');
%   Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
%   EHandle = @STIMA_Lapl_BFE;
%   QuadRule = TProd(gauleg(0,1,NGAUSS));
%   M = assemMat_BFE(Mesh,EHandle,QuadRule);
%   See also Set_Rows, Set_Cols.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);
  
  % Preallocate memory
  
  I = zeros(16*nElements,1);
  J = zeros(16*nElements,1);
  A = zeros(16*nElements,1);
  
  % Assemble element contributions
  
  loc = 1:16;
  for i = 1:nElements
    
    % Extract vertices of current element
    
    idx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(idx,:);
      
    % Compute element contributions
    
    Aloc = EHandle(Vertices,Mesh.ElemFlag(i),varargin{:});
   
    % Add contributions to stiffness matrix
    
    I(loc) = set_Rows(idx,4);
    J(loc) = set_Cols(idx,4);
    A(loc) = Aloc(:);
    loc = loc+16;
    
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