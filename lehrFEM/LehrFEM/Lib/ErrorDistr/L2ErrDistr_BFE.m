function err = L2ErrDistr_BFE(Mesh,u,QuadRule,FHandle,varargin)
% L2ERRDISTR_BFE discretization error distribution in L2 norm for 
% bilinear finite elements.
%
%   ERR = L2ERRDISTR_BFE(MESH,U,QUADRULE,FHANDLE) computes the discretization
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
%   ERR = L2ERRDISTR_BFE(MESH,U,QUADRULE,FHANDLE,FPARAM) also handles the
%   variable length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = L2ErrDistr_BFE(Mesh,u,QuadRule,fhandle);

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  
  nElements = size(Mesh.Elements,1);
    
  % Preallocate memory
  
  err = zeros(nElements,1);
  
  % Precompute shape function values at thhe quedrature points
  
  N = shap_BFE(QuadRule.x);
  grad_N = grad_shap_BFE(QuadRule.x);  

  % Compute discretization error

  for i = 1:nElements
       
    % Extract vertex numbers
    
    vidx = Mesh.Elements(i,:);
    
    % Compute element mapping  
      
    P1 = Mesh.Coordinates(vidx(1),:);
    P2 = Mesh.Coordinates(vidx(2),:);
    P3 = Mesh.Coordinates(vidx(3),:);
    P4 = Mesh.Coordinates(vidx(4),:);
    
    x = N(:,1)*P1+N(:,2)*P2+N(:,3)*P3+N(:,4)*P4;
    z1 = P1(1)*grad_N(:,1:2)+P2(1)*grad_N(:,3:4) + ...
         P3(1)*grad_N(:,5:6)+P4(1)*grad_N(:,7:8);
    z2 = P1(2)*grad_N(:,1:2)+P2(2)*grad_N(:,3:4) + ...
         P3(2)*grad_N(:,5:6)+P4(2)*grad_N(:,7:8);
    det_DPhi_K = abs(z1(:,1).*z2(:,2)-z1(:,2).*z2(:,1));
  
    % Evaluate solutions
      
    u_EX = FHandle(x,varargin{:});
    u_FE = u(vidx(1))*N(:,1)+u(vidx(2))*N(:,2) + ...
           u(vidx(3))*N(:,3)+u(vidx(4))*N(:,4);
      
    % Compute error on current element
      
    err(i) = sum(QuadRule.w.*abs(u_EX-u_FE).^2.*det_DPhi_K);
      
  end
  
  err = sqrt(err);
  
return