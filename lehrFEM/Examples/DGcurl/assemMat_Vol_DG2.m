function varargout = assemMat_Vol_DG2(Mesh,EHandle,varargin)
% ASSEMMAT_VOL_DG2 Assemble vectorial discontinuous Lagrangian nodal FE volume
% contributions.
%
%   A = ASSEMMAT_VOL_DG2(MESH,EHANDLE) assembles the global matrix from the
%   local element contributions given by the function handle EHANDLE and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMAT_VOL_DG2(MESH,EHANDLE) assembles the global matrix
%   from the local element contributions given by the function handle
%   EHANDLE and returns the matrix in an array representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%
%   Example:
%
%   A = assemMat_Vol_DG2(Mesh,@STIMA_Curl);
%  
%   See also set_Rows, set_Cols.

%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  nElements = size(Mesh.Elements,1);  % Number of elements
  
  % Allocate memory
  
  I = zeros(36*nElements,1);
  J = zeros(36*nElements,1);
  A = zeros(36*nElements,1);
  
  % Check for element flags
  
  if(isfield(Mesh,'ElemFlag')),
    ElemFlags = Mesh.ElemFlag; 
  else
    ElemFlags = zeros(nElements,1);
  end
  
  % Assemble element contributions
  
  loc = 1:36;
  for i = 1:nElements
    
    % Extract vertices of current element
    
    Vertices = Mesh.Coordinates(Mesh.Elements(i,:),:);
      
    % Compute element contributions
    
    Aloc = EHandle(Vertices,varargin{:});
   
    % Add contributions to stiffness matrix
    
    idx = 6*(i-1)+[1 2 3 4 5 6];
    
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