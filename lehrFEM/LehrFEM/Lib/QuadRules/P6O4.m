function QuadRule = P6O4()
% P6O4 2D Quadrature rule.
%
%   QUADRULE = P6O4() computes a 6 point quadrature rule of order 4
%   (exact for all polynomials up to degree 3) on the reference element.
%  
%   QUADRULE is a struct containing the following fields:
%    w Weights of the quadrature rule
%    x Abscissae of the quadrature rule (in reference element)
%
%   To recover the barycentric coordinates xbar of the quadrature points
%    xbar = [QuadRule.x, 1-sum(QuadRule.x)'];
%
%   Example:
%
%   QuadRule = P4O3();

%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

 

  QuadRule.w = [1/60; 1/60; 1/60; 9/60; 9/60 ; 9/60];
            
  QuadRule.x = [1/2 0; 1/2 1/2; 0 1/2; 1/6 1/6; 1/6 2/3; 2/3 1/6];
     
return
