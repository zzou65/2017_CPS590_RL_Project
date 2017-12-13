function QuadRule = P7O4()
% P7O4 2D Quadrature rule.
%
%   QUADRULE = P7O4() computes a 7 point quadrature rule of order 4
%   (exact for all polynomials up to degree 3) on the reference element.
%   3/5*P4O3()+2/5*P3O3()
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

 

  QuadRule.w = [1/40; 1/40; 1/40; 2/30; 2/30 ; 2/30;  9/40 ];
            
  QuadRule.x = [0 0; 1 0; 0 1; 1/2 0; 1/2 1/2; 0 1/2; 1/3 1/3];
     
return
