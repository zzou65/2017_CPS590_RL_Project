function err = H1SErr_PDG(Mesh,u,QuadRule,grad_Shap,FHandle,varargin)
% H1SERR_DGCR Discretization error measured in the H1 semi-norm. 
%
%   ERR = H1SERR_DGCR(MESH,U,QUADRULE,FHANDLE) computes the discretization
%   error between the exact solution given by the function handle FHANDLE
%   and the finite element solution U on the struct MESH using the
%   gradients of the reference element shape functions GRAD_SHAP.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh. 
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   GRAD_SHAP is a function pointer to the gradients of the reference
%   element shape functions.
%
%   ERR = H1SERR_DGCR(MESH,U,QUADRULE,GRAD_SHAP,FHANDLE,FPARAM) also
%   handles the variable length argument list FPARAM to the exact solution
%   FHANDLE.
%
%   Example:
%
%   err = H1SErr_PDG(Mesh,u,QuadRule,@grad_shap_DGCR,grad_Uex);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  
  % Check for element flags
  
  if(isfield(Mesh,'ElemFlag')),
    ElemFlag = Mesh.ElemFlag; 
  else
    ElemFlag = zeros(nElements,1);
  end

  % Precompute gradients of shape functions
  
  grad_N = grad_Shap(QuadRule.x);
  nDofs = size(grad_N,2)/2;
  
  % Compute discretization error 
   
  err = 0;
  for i= 1:nElements
        
    % Extract vertex and edge numbers
    
    vidx = Mesh.Elements(i,:);
    idx = nDofs*(i-1)+(1:nDofs);
    
    % Compute element mapping
      
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; ...
          Mesh.Coordinates(vidx(3),:)-bK];
    inv_BK_t = transpose(inv(BK));
    det_BK = abs(det(BK));
    
    % Transform quadrature points
      
    x = QuadRule.x*BK+ones(nPts,1)*bK;
        
    % Evaluate solutions
      
    grad_u_EX = FHandle(x,ElemFlag(i),varargin{:});
    grad_u_FE = zeros(nPts,2);
    for j = 1:nDofs
      loc = 2*(j-1) + [1 2];
      grad_u_FE = grad_u_FE + u(idx(j))*(grad_N(:,loc)*inv_BK_t);
    end
    
    % Compute error on the current element
     
    err = err+sum(QuadRule.w.*sum(abs(grad_u_FE-grad_u_EX).^2,2))*det_BK;
      
  end
  
  err = sqrt(err);
    
return

