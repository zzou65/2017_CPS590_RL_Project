function QuadRule = P1O2()
% P1O2 2D Quadrature rule.
%
%   QUADRULE = P1O2() computes a 1 point Gauss quadrature rule of order 2
%   (exact for all polynomials up to degree 1) on the reference element.
%
%   QUADRULE is a struct containing the following fields:
%    w Weights of the quadrature rule
%    x Abscissae of the quadrature rule
%
%   Example:
%
%   QuadRule = P1O2();

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  QuadRule.w = 1/2;
  QuadRule.x = [1/3 1/3];
  
return