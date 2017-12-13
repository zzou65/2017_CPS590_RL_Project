function varargout = assemMat_LFE2(Mesh,EHandle,varargin)
% ASSEMMAT_LFE2 Assemble nodal FE contributions.
%
%   A = ASSEMMAT_LFE2(MESH,EHANDLE) assembles the global matrix from the
%   local element contributions given by the function handle EHANDLE and
%   returns the matrix in a sparse representation.
%
%   A = ASSEMMAT_LFE2(MESH,EHANDLE,EPARAM) handles the variable length
%   argument list EPARAM to the function handle EHANDLE during the assembly
%   process. 
%
%   [I,J,A] = ASSEMMAT_LFE2(MESH,EHANDLE) assembles the global matrix from
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
%   EHandle = @STIMA_Curl_LFE2;
%   A = assemMat_LFE2(Mesh,EHandle);
%  
%   See also SET_ROWS, SET_COLS.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);
  nCoordinates = size(Mesh.Coordinates,1);
  
  % Preallocate memory
  
  I = zeros(36*nElements,1);
  J = zeros(36*nElements,1);
  A = zeros(36*nElements,1);
  
  % Assemble element contributions
  
  loc = 1:36;
  for i = 1:nElements
    
    % Extract vertices of current element
    
    vidx = Mesh.Elements(i,:);
    idx = [vidx(1) vidx(1)+nCoordinates ...
           vidx(2) vidx(2)+nCoordinates ...
           vidx(3) vidx(3)+nCoordinates];
       
    Vertices = Mesh.Coordinates(vidx,:);
      
    % Compute element contributions
    
    Aloc = EHandle(Vertices,Mesh.ElemFlag(i),varargin{:});
   
    % Add contributions to stiffness matrix
    
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

