function L = assemLoad_Vol_PDG2_vec(Mesh,dim,EHandle,varargin)
%ASSEMLOAD_VOL_PDG2_VEC Assemble vector-valued DG volume contributions
%
%   L = ASSEM_LOAD_VOL_PDG2_VEC(MESH,DIM,EHANDLE) assembles the global load
%   vector from the local element contributions given by the function
%   handle EHANDLE.  This code is for vector-valued functions, for example
%   mixed DG, with the same basis functions in each component of the image.
%   The structure of the discrete vectors should be such that all
%   coefficients corresponding to the first dimension come first, then the 
%   second dimension, and so on.
%
%   L = ASSEM_LOAD_VOL_PDG2_VEC(MESH,DIM,EHANDLE,EPARAM) passes the
%   variable-length argument list EPARAM to the function handle EHANDLE
%   during the assembly process.
%
%   The struct MESH must contain at least the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMDATA     N-by-1 structure array containing at least the fields:
%       NDOFS       The number of degrees of freedom on the corresponding
%                   element.

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
  nDofsTot = nDofsSum(end);         % Total number of degrees of freedom per dimension
  numel = dim*nDofsTot;             % Number of matrix entries
  
  % Preallocate memory
  L = zeros(numel,1);
  
  % Assemble element contributions
  dimOffset = nDofsTot*(0:dim-1);
  for i = 1:nElements
    
    % Extract vertices
    vidx = Mesh.Elements(i,:);
    
    % Compute element contributions
    Lloc = EHandle(Mesh.Coordinates(vidx,:),Mesh.ElemData(i),varargin{:});
    
    % Add contributions to global load vector
    idx0 = nDofsSum(i)+(1:nDofs(i));
    idx = zeros(1,dim*nDofs(i));
    idx(:) = idx0(ones(dim,1),:)' + dimOffset(ones(nDofs(i),1),:);
    L(idx) = L(idx)+Lloc;
    
  end
  
return