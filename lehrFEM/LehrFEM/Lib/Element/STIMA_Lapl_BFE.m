function Aloc = STIMA_Lapl_BFE(Vertices,ElemInfo,QuadRule,varargin)
% STIMA_LAPL_BFE Element stiffness matrix for the Laplacian.
%
%   ALOC = STIMA_LAPL_BFE(VERTICES) computes the element stiffness matrix
%   for the Laplacian using bi-linear Lagrangian finite elements.
%
%   VERTICES is 4-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   Example:
%
%   Aloc=STIMA_Lapl_BFE([0 0; 1 0; 1 1;0 1],0,gauleg_sqr(5));

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nGauss = size(QuadRule.w,1);
  
  % Preallocate memory
  
  Aloc = zeros(4,4);
  DPhi_K = zeros(2,2);
  
  % Compute values of shape functions
  
  grad_N = grad_shap_BFE(QuadRule.x);
  
  % Extract vertices
  
  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  P4 = Vertices(4,:);
  
  for i = 1:nGauss

    % Compute element map  
    
    DPhi_K(1,:) = P1(1)*grad_N(i,1:2)+P2(1)*grad_N(i,3:4) + ...
                  P3(1)*grad_N(i,5:6)+P4(1)*grad_N(i,7:8);
    DPhi_K(2,:) = P1(2)*grad_N(i,1:2)+P2(2)*grad_N(i,3:4) + ...
                  P3(2)*grad_N(i,5:6)+P4(2)*grad_N(i,7:8);
    inv_DPhi_K = inv(DPhi_K);
    TK =inv_DPhi_K*transpose(inv_DPhi_K)*abs(det(DPhi_K));
        
    % Compute entries of stiffness matrix
    
    Aloc(1,1) = Aloc(1,1) + QuadRule.w(i)*grad_N(i,1:2)*TK*transpose(grad_N(i,1:2));
    Aloc(1,2) = Aloc(1,2) + QuadRule.w(i)*grad_N(i,1:2)*TK*transpose(grad_N(i,3:4));
    Aloc(1,3) = Aloc(1,3) + QuadRule.w(i)*grad_N(i,1:2)*TK*transpose(grad_N(i,5:6));
    Aloc(1,4) = Aloc(1,4) + QuadRule.w(i)*grad_N(i,1:2)*TK*transpose(grad_N(i,7:8));
    Aloc(2,2) = Aloc(2,2) + QuadRule.w(i)*grad_N(i,3:4)*TK*transpose(grad_N(i,3:4));
    Aloc(2,3) = Aloc(2,3) + QuadRule.w(i)*grad_N(i,3:4)*TK*transpose(grad_N(i,5:6));
    Aloc(2,4) = Aloc(2,4) + QuadRule.w(i)*grad_N(i,3:4)*TK*transpose(grad_N(i,7:8));
    Aloc(3,3) = Aloc(3,3) + QuadRule.w(i)*grad_N(i,5:6)*TK*transpose(grad_N(i,5:6));
    Aloc(3,4) = Aloc(3,4) + QuadRule.w(i)*grad_N(i,5:6)*TK*transpose(grad_N(i,7:8));
    Aloc(4,4) = Aloc(4,4) + QuadRule.w(i)*grad_N(i,7:8)*TK*transpose(grad_N(i,7:8));
    
  end
  
  % Fill in lower triangular part
  
  tri = triu(Aloc);
  Aloc = tri+tril(tri',-1);
  
return