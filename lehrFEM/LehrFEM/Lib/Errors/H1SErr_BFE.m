function err = H1SErr_BFE(Mesh,u,QuadRule,FHandle,varargin)
% H1SERR_BFE Discretization error in H1 semi-norm for bi-linear finite elements.
%
%   ERR = HISERR_BFE(MESH,U,QUADRULE,FHANDLE) computes the discretization
%   error between the exact solution given by the function handle FHANDLE
%   and the finite element solution U on the struct MESH.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-4 matrix specifying the elements of the mesh. 
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   ERR = HISERR_BFE(MESH,U,QUADRULE,FHANDLE,FPARAM) also handles the
%   variable length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = HISErr_BFE(Mesh,u,QuadRule,FHandle);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);
  nElements = size(Mesh.Elements,1);
  
  % Precompute shap functions and gradients of shape functions
  
  grad_N = grad_shap_BFE(QuadRule.x);
  N = shap_BFE(QuadRule.x);
  
  % Compute discretization error 
   
  err = 0;
  DPhi_K = zeros(2,2);
  for i= 1:nElements
        
    % Extract vertex numbers
    
    vidx = Mesh.Elements(i,:);
         
    % Compute element mapping  
      
    P1 = Mesh.Coordinates(vidx(1),:);
    P2 = Mesh.Coordinates(vidx(2),:);
    P3 = Mesh.Coordinates(vidx(3),:);
    P4 = Mesh.Coordinates(vidx(4),:);
    
    % Evaluate solutions
    
    x = N(:,1)*P1+N(:,2)*P2+N(:,3)*P3+N(:,4)*P4;
    grad_u_EX = FHandle(x,varargin{:});
    
    for j = 1:nPts
      DPhi_K(1,:) = P1(1)*grad_N(j,1:2)+P2(1)*grad_N(j,3:4) + ...
                    P3(1)*grad_N(j,5:6)+P4(1)*grad_N(j,7:8);
      DPhi_K(2,:) = P1(2)*grad_N(j,1:2)+P2(2)*grad_N(j,3:4) + ...
                    P3(2)*grad_N(j,5:6)+P4(2)*grad_N(j,7:8);
      inv_DPhi_K = inv(DPhi_K);
      det_DPhi_K = abs(det(DPhi_K));       
        
      grad_u_FE = (u(vidx(1))*grad_N(j,1:2)+ ...
                   u(vidx(2))*grad_N(j,3:4)+ ...
                   u(vidx(3))*grad_N(j,5:6)+ ...
                   u(vidx(4))*grad_N(j,7:8))*inv_DPhi_K;
    
      err = err + QuadRule.w(j)*sum((grad_u_EX(j,:)-grad_u_FE).^2)*det_DPhi_K;
 
    end   
  end
  err = sqrt(err);
    
return