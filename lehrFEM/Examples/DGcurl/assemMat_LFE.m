function varargout = assemMat_LFE(Mesh,EHandle,varargin)
% ASSEMMAT_LFE Assemble linear FE contributions.
%
%   A = ASSEMMAT_LFE(MESH,EHANDLE) assembles the global matrix from the
%   local element contributions given by the function handle EHANDLE and
%   returns the matrix in a sparse representation.
%
%   A = ASSEMMAT_LFE(MESH,EHANDLE,EPARAM) handles the variable length
%   argument list EPARAM to the function handle EHANDLE during the assembly
%   process. 
%
%   [I,J,A] = ASSEMMAT_LFE(MESH,EHANDLE) assembles the global matrix from
%   the local element contributions given by the function handle EHANDLE
%   and returns the matrix in an array representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 or N-by-4 matrix specifying the elements of the
%                 mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
%   EHandle = @STIMA_Lapl_LFE;
%   A = assemMat_LFE(Mesh,EHandle);
%  
%   See also set_Rows, set_Cols.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);
  
  % Preallocate memory
  
  I = zeros(9*nElements,1);
  J = zeros(9*nElements,1);
  A = zeros(9*nElements,1);
  
  % Check for element flags
  
  if(isfield(Mesh,'ElemFlag')),
    flags = Mesh.ElemFlag; 
  else
    flags = zeros(nElements,1);
  end

  % Assemble element contributions
  
  loc = 1:9;
  for i = 1:nElements
    
    % Extract vertices of current element
    
    idx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(idx,:);
      
    % Compute element contributions
    
    Aloc = EHandle(Vertices,flags(i),varargin{:});
   
    % Add contributions to stiffness matrix
    
    I(loc) = set_Rows(idx,3);
    J(loc) = set_Cols(idx,3);
    A(loc) = Aloc(:);
    loc = loc+9;
    
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