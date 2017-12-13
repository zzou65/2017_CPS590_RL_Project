function QuadRule = P3O3()
% P3O3 2D Quadrature rule.
%
%   QUADRULE = P3O3() computes a 3 point Gauss quadrature rule of order 3
%   (exact for all polynomials up to degree 2) on the reference element.
%
%   QUADRULE is a struct containing the following fields:
%    w Weights of the quadrature rule
%    x Abscissae of the quadrature rule
%
%   Example:
%
%   QuadRule = P3O3();
%

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  QuadRule.w = 1/6*ones(3,1);
  QuadRule.x = [  0 1/2; ...
                1/2 1/2;
                1/2   0];
  
return