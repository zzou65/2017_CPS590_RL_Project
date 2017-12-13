function varargout = assemMat_Vol_PDG2_vec(Mesh,dim,EHandle,varargin)
%ASSEMMAT_VOL_PDG2_VEC Assemble vector-valued DG volume contributions
%
%   A = ASSEMMAT_VOL_PDG2_VEC(MESH,DIM,EHANDLE) assembles the global matrix
%   from the local element contributions given by the function handle
%   EHANDLE and returns the matrix in sparse representation.  This code is
%   for vector-valued functions, for example mixed DG, with the same basis
%   functions in each component of the image.  The structure of the
%   discrete vectors should be such that all coefficients corresponding to
%   the first dimension come first, then the second dimension, and so on.
%
%   A = ASSEMMAT_VOL_PDG2_VEC(MESH,DIM,EHANDLE,EPARAM) passes the variable-
%   length argument list EPARAM to the function handle EHANDLE during the
%   assembly process.
%
%   [I,J,A] = ASSEMMAT_VOL_PDG2_VEC(...) returns the matrix in an array
%   representation.
%
%   DIM is the dimension of the image space.  For example, for mixed DG,
%   the image space consists of one dimension for the function value and
%   two for the gradient, amounting to a total of three.  If the matrix is
%   not square, use DIM = [DIM1,DIM2], where DIM1 is the dimension of the
%   elements of the space determining the number of rows of the matrix (ie.
%   the dimension of the test functions in finite elements) and DIM2 is the
%   dimension of the elements of the space determining the number of
%   columns of the matrix (ie. the dimension of the solution).
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

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  nElements = size(Mesh.Elements,1); % Number of elements
  dim1 = dim(1);
  dim2 = dim(end);

  % Count degrees of freedom on all elements and calculate number of
  % matrix entries
  nDofs = [Mesh.ElemData.nDofs];    % Number of degrees of freedom on element
  nDofsSum = cumsum([0,nDofs]);     % Number of degrees of freedom on all preceding elements
  nDofsTot = nDofsSum(end);         % Total number of degrees of freedom per dimension
  numel = dim1*dim2*sum(nDofs.^2);  % Number of matrix entries
  
  % Preallocate memory
  I = zeros(numel,1);
  J = zeros(numel,1);
  A = zeros(numel,1);
  
  % Assemble element contributions
  loc = 1;
  dim1Offset = nDofsTot*(0:dim1-1);
  dim2Offset = nDofsTot*(0:dim2-1);
  for i=1:nElements
    
    % Extract vertices of current element
    Vertices = Mesh.Coordinates(Mesh.Elements(i,:),:);
      
    % Compute element contributions
    Aloc = EHandle(Vertices,Mesh.ElemData(i),varargin{:});
   
    % Add contributions to stiffness matrix
    idx0 = nDofsSum(i)+(1:nDofs(i));
    idx1 = zeros(1,dim1*nDofs(i));
    idx1(:) = idx0(ones(dim1,1),:)' + dim1Offset(ones(nDofs(i),1),:);
    idx2 = zeros(1,dim2*nDofs(i));
    idx2(:) = idx0(ones(dim2,1),:)' + dim2Offset(ones(nDofs(i),1),:);
    loc1 = loc + dim1*dim2*nDofs(i)^2 - 1;
   
    I(loc:loc1) = set_Rows(idx1,dim2*nDofs(i));
    J(loc:loc1) = set_Cols(idx2,dim1*nDofs(i));
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