function Mloc = MASS_BFE(Vertices,ElemInfo,QuadRule,varargin)
% MASS_BFE Element mass matrix.
%
%   M = MASS_BFE(VERTICES,ELEMINFO,QUADRULE) computes the element mass 
%   matrix using bi-linear Lagrangian finite elements.
%
%   VERTICES is a 4-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   ELEMINFO is an integer parameter which is used to specify additional
%   element information on each element.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used to
%   do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%   
%
%   Example:
%
%   M = MASS_BFE([0 0; 1 0; 1 1; 0 1],0,TProd(gauleg(0,1,NGAUSS)));
%
%   See also LOAD_BFE.

%   Copyright 2005-2005 Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory
  
  Mloc = zeros(4,4);
  
  % Compute element mapping
  
  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  P4 = Vertices(4,:);
  
  N = shap_BFE(QuadRule.x);
  grad_N = grad_shap_BFE(QuadRule.x);
  
  z1 = P1(1)*grad_N(:,1:2)+P2(1)*grad_N(:,3:4) + ...
       P3(1)*grad_N(:,5:6)+P4(1)*grad_N(:,7:8);
  z2 = P1(2)*grad_N(:,1:2)+P2(2)*grad_N(:,3:4) + ...
       P3(2)*grad_N(:,5:6)+P4(2)*grad_N(:,7:8);
  det_DPhi_K = abs(z1(:,1).*z2(:,2)-z1(:,2).*z2(:,1));
  
  % Compute the matrix
  
  for i = 1:4
    for j = i:4
      Mloc(i,j) = sum(QuadRule.w.*N(:,i).*N(:,j).*det_DPhi_K); 
    end
  end
  
  % Fill in lower triangular part
  
  tri = triu(Mloc);
  Mloc = tri+tril(tri',-1);
  
return