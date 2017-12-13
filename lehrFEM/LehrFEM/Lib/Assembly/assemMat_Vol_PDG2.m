function varargout = assemMat_Vol_PDG2(Mesh,EHandle,varargin)
%ASSEMMAT_VOL_PDG2 Assemble DG volume contributions
%
%   A = ASSEMMAT_VOL_PDG2(MESH,EHANDLE) assembles the global matrix from
%   the local element contributions given by the function handle EHANDLE
%   and returns the matrix in sparse representation.
%
%   A = ASSEMMAT_VOL_PDG2(MESH,EHANDLE,EPARAM) passes the variable-length
%   argument list EPARAM to the function handle EHANDLE during the assembly
%   process.
%
%   [I,J,A] = ASSEMMAT_VOL_PDG2(...) returns the matrix in an array
%   representation.
%
%   EHANDLE must take the arguments Vertices and ElemData, where Vertices
%   is a matrix containing the coordinates of the vertices of the element
%   and ElemData is a struct containing data relevant to the element, such
%   as the number of degrees of freedom on the element.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 or N-by-4 matrix specifying the elements of the
%                 mesh.
%    ELEMDATA     N-by-1 structure array containing at least the fields:
%       NDOFS       The number of degrees of freedom on the corresponding
%                   element.
%
%   See also assemMat_Vol_PDG, assemMat_Inn_PDG2, assemMat_Bnd_PDG2.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  nElements = size(Mesh.Elements,1); % Number of elements

  % Count degrees of freedom on all elements and calculate number of
  % matrix entries
  nDofs = [Mesh.ElemData.nDofs];    % Number of degrees of freedom on element
  nDofsSum = cumsum([0,nDofs]);     % Number of degrees of freedom on all preceding elements
  numel = sum(nDofs.^2);            % Number of matrix entries
  
  % Preallocate memory
  I = zeros(numel,1);
  J = zeros(numel,1);
  A = zeros(numel,1);
  
  % Assemble element contributions
  loc = 1;
  for i=1:nElements
    
    % Extract vertices of current element
    Vertices = Mesh.Coordinates(Mesh.Elements(i,:),:);
      
    % Compute element contributions
    Aloc = EHandle(Vertices,Mesh.ElemData(i),varargin{:});
   
    % Add contributions to stiffness matrix
    idx = nDofsSum(i)+(1:nDofs(i));
    loc1 = loc + nDofs(i)^2 - 1;
   
    I(loc:loc1) = set_Rows(idx,nDofs(i));
    J(loc:loc1) = set_Cols(idx,nDofs(i));
    A(loc:loc1) = Aloc(:);
    
    loc = loc1 + 1;
      
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