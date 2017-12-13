function varargout = assemMat_P1P0(Mesh,EHandle,varargin)
% ASSEMMAT_P1P0 Assemble linear/constant FE contributions.
%
%   A = ASSEMMAT_P1P0(MESH,EHANDLE) assembles the global matrix from the
%   local element contributions given by the function handle EHANDLE and
%   returns the matrix in a sparse representation.
%
%   A = ASSEMMAT_P1P0(MESH,EHANDLE,EPARAM) handles the variable length
%   argument list EPARAM to the function handle EHANDLE during the assembly
%   process. 
%
%   [I,J,A] = ASSEMMAT_P1P0(MESH,EHANDLE) assembles the global matrix from
%   the local element contributions given by the function handle EHANDLE
%   and returns the matrix in an array representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
%   EHandle = @MASS_P1P0;
%   A = assemMat_P1P0(Mesh,EHandle);
%  
%   See also set_Rows, set_Cols.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);
  
  % Preallocate memory
  
  I = zeros(3*nElements,1);
  J = zeros(3*nElements,1);
  A = zeros(3*nElements,1);
  
  % Assemble element contributions
  
  loc = 1:3;
  for i = 1:nElements
    
    % Extract vertices of current element
    
    vidx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(vidx,:);
      
    % Compute element contributions
    
    Aloc = EHandle(Vertices,Mesh.ElemFlag(i),varargin{:});
   
    % Add contributions to stiffness matrix
    
    row_idx = vidx;
    col_idx = i;
    I(loc) = set_Rows(row_idx,1);
    J(loc) = set_Cols(col_idx,3);
    A(loc) = Aloc(:);
    loc = loc+3;
    
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