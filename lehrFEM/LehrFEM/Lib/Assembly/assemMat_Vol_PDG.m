function varargout = assemMat_Vol_PDG(Mesh,nDofs,EHandle,varargin)
% ASSEMMAT_VOL_PDG Assemble discontinuous volume contributions.
%
%   A = ASSEMMAT_VOL_PDG(MESH,NDOFS,EHANDLE) assembles the global matrix
%   from the local element contributions given by the function handle
%   EHANDLE and returns the matrix in a sparse representation. The integer
%   NDOFS specifies the number of dofs per element.
%
%   A = ASSEMMAT_VOL_PDG(MESH,NDOFS,EHANDLE,EPARAM) handles the variable
%   length argument list EPARAM to the function handle EHANDLE during the
%   assembly process. 
%
%   [I,J,A] = ASSEMMAT_VOL_DG(MESH,NDOFS,EHANDLE) assembles the global
%   matrix from the local element contributions given by the function
%   handle EHANDLE and returns the matrix in an array representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 or N-by-4 matrix specifying the elements of the
%                 mesh.
%
%   Example:
%
%   A = assemMat_Vol_PDG(Mesh,3,@MASS_Vol_PDG,P3O3(),@shap_DGCR);
%  
%   See also set_Rows, set_Cols.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);

  % Preallocate memory
  
  I = zeros(nDofs^2*nElements,1);
  J = zeros(nDofs^2*nElements,1);
  A = zeros(nDofs^2*nElements,1);
  
  % Check for element flags
  
  if(isfield(Mesh,'ElemFlag')),
    ElemFlags = Mesh.ElemFlag; 
  else
    ElemFlags = zeros(nElements,1);
  end
  
  % Assemble element contributions
  
  loc = 1:(nDofs^2);
  for i = 1:nElements
    
    % Extract vertices of current element
    
    Vertices = Mesh.Coordinates(Mesh.Elements(i,:),:);
      
    % Compute element contributions
    
    Aloc = EHandle(Vertices,ElemFlags(i),varargin{:});
   
    % Add contributions to stiffness matrix
    
    idx = nDofs*(i-1)+(1:nDofs);
   
    I(loc) = set_Rows(idx,nDofs);
    J(loc) = set_Cols(idx,nDofs);
    A(loc) = Aloc(:);
    
    loc = loc+nDofs^2;
      
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
