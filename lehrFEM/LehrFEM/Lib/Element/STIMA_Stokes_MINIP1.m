function Aloc = STIMA_Stokes_MINIP1(Vertices,ElemInfo,nu,QuadRule,varargin)
% STIMA_STOKES_MINIP1 Element stiffness matrix for the stokes problem.
%
%   ALOC = STIMA_STOKES_MINIP1(VERTICES,ELEMINFO,NU,QUADRULE) computes the
%   element stiffness matrix for the steady Stokes problem with the
%   viscosity NU.
%
%   VERTICES is a 3-by-2 matrix specifying the vertices of the current
%   element in a row wise orientation.
%
%   ELEMINFO is an integer parameter which is used to specify additional
%   element information on each element.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%
%   Example:
%
%   nu = 1;   
%   Aloc = STIMA_Stokes_MINIP1([0 0; 1 0; 0 1],0,nu,P7O6());
%
%   See also grad_shap_MINI, shap_LFE.

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory
  
  Aloc = zeros(12,12);
  
  % Compute element mapping
  
  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK; Vertices(3,:)-bK];
  inv_BK_t = transpose(inv(BK));
  det_BK = abs(det(BK));
   
  % Compute element shape functions
   
  grad_N = grad_shap_MINI(QuadRule.x);
  grad_N(:,1:2) = grad_N(:,1:2)*inv_BK_t;
  grad_N(:,3:4) = grad_N(:,3:4)*inv_BK_t;
  grad_N(:,5:6) = grad_N(:,5:6)*inv_BK_t;
  grad_N(:,7:8) = grad_N(:,7:8)*inv_BK_t;
  N = shap_LFE(QuadRule.x);
  
  % Compute stiffness matrix
  
  Aloc(1,1) = nu*sum(QuadRule.w.*sum(grad_N(:,1:2).*grad_N(:,1:2),2))*det_BK;
  Aloc(1,2) = nu*sum(QuadRule.w.*sum(grad_N(:,1:2).*grad_N(:,3:4),2))*det_BK;
  Aloc(1,3) = nu*sum(QuadRule.w.*sum(grad_N(:,1:2).*grad_N(:,5:6),2))*det_BK;
  Aloc(1,4) = nu*sum(QuadRule.w.*sum(grad_N(:,1:2).*grad_N(:,7:8),2))*det_BK;
  Aloc(2,2) = nu*sum(QuadRule.w.*sum(grad_N(:,3:4).*grad_N(:,3:4),2))*det_BK;
  Aloc(2,3) = nu*sum(QuadRule.w.*sum(grad_N(:,3:4).*grad_N(:,5:6),2))*det_BK;
  Aloc(2,4) = nu*sum(QuadRule.w.*sum(grad_N(:,3:4).*grad_N(:,7:8),2))*det_BK;
  Aloc(3,3) = nu*sum(QuadRule.w.*sum(grad_N(:,5:6).*grad_N(:,5:6),2))*det_BK;
  Aloc(3,4) = nu*sum(QuadRule.w.*sum(grad_N(:,5:6).*grad_N(:,7:8),2))*det_BK;
  Aloc(4,4) = nu*sum(QuadRule.w.*sum(grad_N(:,7:8).*grad_N(:,7:8),2))*det_BK;
  
  Aloc(1,9) = sum(QuadRule.w.*grad_N(:,1).*N(:,1))*det_BK;
  Aloc(1,10) = sum(QuadRule.w.*grad_N(:,1).*N(:,2))*det_BK;
  Aloc(1,11) = sum(QuadRule.w.*grad_N(:,1).*N(:,3))*det_BK;
  Aloc(2,9) = sum(QuadRule.w.*grad_N(:,3).*N(:,1))*det_BK;
  Aloc(2,10) = sum(QuadRule.w.*grad_N(:,3).*N(:,2))*det_BK;
  Aloc(2,11) = sum(QuadRule.w.*grad_N(:,3).*N(:,3))*det_BK;
  Aloc(3,9) = sum(QuadRule.w.*grad_N(:,5).*N(:,1))*det_BK;
  Aloc(3,10) = sum(QuadRule.w.*grad_N(:,5).*N(:,2))*det_BK;
  Aloc(3,11) = sum(QuadRule.w.*grad_N(:,5).*N(:,3))*det_BK;
  Aloc(4,9) = sum(QuadRule.w.*grad_N(:,7).*N(:,1))*det_BK;
  Aloc(4,10) = sum(QuadRule.w.*grad_N(:,7).*N(:,2))*det_BK;
  Aloc(4,11) = sum(QuadRule.w.*grad_N(:,7).*N(:,3))*det_BK;
  
  Aloc(5,5) = Aloc(1,1);
  Aloc(5,6) = Aloc(1,2);
  Aloc(5,7) = Aloc(1,3);
  Aloc(5,8) = Aloc(1,4);
  Aloc(6,6) = Aloc(2,2);
  Aloc(6,7) = Aloc(2,3);
  Aloc(6,8) = Aloc(2,4);
  Aloc(7,7) = Aloc(3,3);
  Aloc(7,8) = Aloc(3,4);
  Aloc(8,8) = Aloc(4,4);
  
  Aloc(5,9) = sum(QuadRule.w.*grad_N(:,2).*N(:,1))*det_BK;
  Aloc(5,10) = sum(QuadRule.w.*grad_N(:,2).*N(:,2))*det_BK;
  Aloc(5,11) = sum(QuadRule.w.*grad_N(:,2).*N(:,3))*det_BK;
  Aloc(6,9) = sum(QuadRule.w.*grad_N(:,4).*N(:,1))*det_BK;
  Aloc(6,10) = sum(QuadRule.w.*grad_N(:,4).*N(:,2))*det_BK;
  Aloc(6,11) = sum(QuadRule.w.*grad_N(:,4).*N(:,3))*det_BK;
  Aloc(7,9) = sum(QuadRule.w.*grad_N(:,6).*N(:,1))*det_BK;
  Aloc(7,10) = sum(QuadRule.w.*grad_N(:,6).*N(:,2))*det_BK;
  Aloc(7,11) = sum(QuadRule.w.*grad_N(:,6).*N(:,3))*det_BK;
  Aloc(8,9) = sum(QuadRule.w.*grad_N(:,8).*N(:,1))*det_BK;
  Aloc(8,10) = sum(QuadRule.w.*grad_N(:,8).*N(:,2))*det_BK;
  Aloc(8,11) = sum(QuadRule.w.*grad_N(:,8).*N(:,3))*det_BK;
  
  Aloc(9,12) = sum(QuadRule.w.*N(:,1))*det_BK;
  Aloc(10,12) = sum(QuadRule.w.*N(:,2))*det_BK;
  Aloc(11,12) = sum(QuadRule.w.*N(:,3))*det_BK;
  
  % Fill in lower triangular part
  
  tri = triu(Aloc);
  Aloc = tri+tril(tri',-1);
  
return