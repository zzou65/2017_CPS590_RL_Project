function [PL,PR] = assemPrec_SVD_PDG(Mesh,B)
%ASSEMPREC_SVD_PDG Construct preconditioner for arbitrary DG method
%
%   [PL,PR] = ASSEMPREC_SVD_PDG(NDOFS,B) constructs left and right
%   preconditioners PL and PR based on SVD of the block diagonal of B with
%   block sizes given by NDOFS.  The blocks represent degrees of freedom on
%   the same element.  If the number of degrees of freedom is constant,
%   NDOFS can be a scalar; otherwise, it should be a vector containing the
%   number of degrees of freedom on each element.
%
%   [PL,PR] = ASSEMPREC_SVD_PDG(MESH,B) constructs left and right
%   preconditioners PL and PR based on SVD of the block diagonal of B with
%   block sizes given by MESH.ElemData(:).nDofs.  The struct MESH must 
%   contain at least the following fields:
%    ELEMENTS     N-by-k matrix listing the elements of the mesh
%    ELEMDATA     N-by-1 structure array containing at least the fields:
%       NDOFS       The number of degrees of freedom on the corresponding
%                   element.
%
%   The matrix B can be a stiffness matrix or mass matrix or any other
%   reasonable matrix of the correct dimension.  Its block diagonal must,
%   however, be regular.
%
%   The matrices PL and PR are sparse block diagonal such that for any 
%   matrix A = B+E for small E,
%     PL*A*PR
%   is well-conditioned.  An equation A*x=b can be solved in the following
%   two steps:
%     y = (PL*A*PR)\(PL*b);
%     x = PR*y;
%
%   More precisely, for each block B_loc of B, let [U,S,V] be the singular
%   value decomposition of B_loc.  If D is the inverse of the square root
%   of S, then the corresponding blocks of the preconditioners are
%   PL_loc=(U*D)' and PR_loc=V*D.  In particular, PL_loc*B_loc*PR_loc is
%   the unit matrix.
%
%   See also svd.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize Constants
  if(isstruct(Mesh))
    nElements = size(Mesh.Elements,1); % Number of elements in mesh
    nDofs = [Mesh.ElemData.nDofs];     % Number of degrees of freedom on element
  elseif(isscalar(Mesh))
    nElements = floor(size(B,1)/Mesh); % Number of elements in mesh
    nDofs = Mesh(ones(1,nElements));   % Number of degrees of freedom on element
  elseif(isnumeric(Mesh))
    nElements = length(Mesh);          % Number of elements in mesh
    nDofs = Mesh;                      % Number of degrees of freedom on element
  end
  nDofsSum = cumsum([0,nDofs]);      % Number of degrees of freedom on all preceding elements
  numel = sum(nDofs.^2);             % Number of matrix entries
  
  % Prealocate memory
  I = zeros(numel,1);
  J = zeros(numel,1);
  Pl = zeros(numel,1);
  Pr = zeros(numel,1);
  
  % Loop over elements
  loc = 1;
  for j=1:nElements
    
    % Compute indices
    idx = nDofsSum(j)+(1:nDofs(j));
    loc1 = loc + nDofs(j)^2 - 1;
    
    % Compute SVD of block of B
    Bloc = full(B(idx,idx));
    [U,S,V] = svd(Bloc);
    
    % Compute local basis transformations in domain and codomain
    D = diag(sqrt(1./diag(S)));
    U = (U*D)';
    V = V*D;
    
    % Add contributions to global matrices
    I(loc:loc1) = set_Rows(idx,nDofs(j));
    J(loc:loc1) = set_Cols(idx,nDofs(j));
    Pl(loc:loc1) = U(:);
    Pr(loc:loc1) = V(:);
    
    % Update index
    loc = loc1 + 1;
    
  end
  
  % Assemble preconditioner matrices
  PL = sparse(I,J,Pl);
  PR = sparse(I,J,Pr);

return