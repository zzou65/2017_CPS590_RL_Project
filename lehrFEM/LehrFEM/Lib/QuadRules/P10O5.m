function QuadRule = P10O5()
% P7O4 2D Quadrature rule.
%
%   QUADRULE = P10O5() computes a 10 point quadrature rule of order 5
%   (exact for all polynomials up to degree 4) on the reference element.
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

 

  QuadRule.w = [1/120; 1/120; 1/120; 1/30; 1/30 ; 1/30;...
	1/10; 1/10; 1/10; 3/40];
            
  QuadRule.x = [0 0; 1 0; 0 1;...
	1/2 0; 1/2 1/2; 0 1/2;...
        1/6 1/6; 1/6 2/3; 2/3 1/6;...	
	1/3 1/3];
     
return
