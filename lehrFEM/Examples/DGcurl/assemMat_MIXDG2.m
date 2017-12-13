function varargout = assemMat_MIXDG2(Mesh,EHandle,varargin)
% ASSEMMAT_MIXDG2 Assemble nodal FE contributions.
%
%   A = ASSEMMAT_MIXDG2(MESH,EHANDLE) assembles the global matrix from the
%   local element contributions given by the function handle EHANDLE and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMAT_MIXDG2(MESH,EHANDLE) assembles the global matrix from
%   the local element contributions given by the function handle EHANDLE
%   and returns the matrix in an array representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 or N-by-4 matrix specifying the elements of the
%                 mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.

%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);
  nCoordinates = size(Mesh.Coordinates,1);
  
  % Preallocate memory
  
  I = zeros(18*nElements,1);
  J = zeros(18*nElements,1);
  A = zeros(18*nElements,1);
  
  if(isfield(Mesh,'ElemFlag')),
    ElemFlags = Mesh.ElemFlag; 
  else
    ElemFlags = zeros(nElements,1);
  end
  
  % Assemble element contributions
  
  loc = 1:18;
  for i = 1:nElements
    
    % Extract vertices of current element
    
    vidx = Mesh.Elements(i,:);
    idx = [vidx(1) ...
           vidx(2) ...
           vidx(3) ];
       
    Vertices = Mesh.Coordinates(vidx,:);
      
    % Compute element contributions
    
    Aloc = EHandle(Vertices,Mesh.ElemFlag(i),varargin{:});
    
    % Add contributions to stiffness matrix
    
    idx_j = 3*(i-1)+[1 2 3];
    idx_i = 6*(i-1)+[1 2 3 4 5 6];
    J(loc) = set_Cols(idx_j,6);
    I(loc) = set_Rows(idx_i,3);
    A(loc) = Aloc(:);
    loc = loc+18;
    
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

