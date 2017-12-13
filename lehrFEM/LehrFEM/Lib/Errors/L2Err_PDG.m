function err = L2Err_PDG(Mesh,u,QuadRule,Shap,FHandle,varargin)
% L2ERR_PDG Discretization error measured in the L2 norm.
%
%   ERR = L2ERR_PDG(MESH,U,QUADRULE,SHAP,FHANDLE) computes the
%   discretization error between the exact solution given by the function
%   handle FHANDLE and the finite element solution U on the struct MESH
%   using the reference element shape functions SHAP.
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
%   SHAP is a function pointer to the reference element shape functions.
%
%   ERR = L2ERR_PDG(MESH,U,QUADRULE,SHAP,FHANDLE,FPARAM) also handles the
%   variable length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = L2Err_PDG(Mesh,u,QuadRule,@shap_DGCR,Uex);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  
  nPts = size(QuadRule.w,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  
  % Check for element flags
  
  if(isfield(Mesh,'ElemFlag')),
    ElemFlag = Mesh.ElemFlag; 
  else
    ElemFlag = zeros(nElements,1);
  end

  % Precompute shape function values at the quedrature points
  
  N = Shap(QuadRule.x);
  nDofs = size(N,2);  
  
  % Compute discretization error

  err = 0;
  for i = 1:nElements
       
    % Extract dof numbers on the current element
    
    idx = nDofs*(i-1)+(1:nDofs);
          
    % Compute element mapping  
        
    bK = Mesh.Coordinates(Mesh.Elements(i,1),:);
    BK = [Mesh.Coordinates(Mesh.Elements(i,2),:)-bK; ...
          Mesh.Coordinates(Mesh.Elements(i,3),:)-bK];
    det_BK = abs(det(BK));

    % Transform quadrature points
      
    x = QuadRule.x*BK+ones(nPts,1)*bK;
      
    % Evaluate solutions
      
    u_EX = FHandle(x,ElemFlag(i),varargin{:});
    u_FE = zeros(nPts,1);
    for j = 1:nDofs
      u_FE = u_FE + u(idx(j))*N(:,j);  
    end
    
    % Compute error on current element
      
    err = err + sum(QuadRule.w.*abs(u_EX-u_FE).^2)*det_BK;
      
  end
  
  err = sqrt(err);
  
return