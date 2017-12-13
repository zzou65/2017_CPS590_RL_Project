 function QuadRule = P10O4()
 % P7O4 2D Quadrature rule.
 %
 %   QUADRULE = P10O5() computes a 10 point quadrature rule of order

 %   (exact for all polynomials up to degree 4) on the reference

 %
 %   QUADRULE is a struct containing the following fields:
 %    w Weights of the quadrature rule
 %    x Abscissae of the quadrature rule (in reference

 %
 %   To recover the barycentric coordinates xbar of the

 %    xbar = [QuadRule.x, 1-sum(QuadRule.x)'];
 %
 %   Example:
 %
 %   QuadRule = P4O3();

 %   Copyright 2008-2008 Holger Heumann
 %   SAM - Seminar for Applied Mathematics
 %   ETH-Zentrum
 %   CH-8092 Zurich, Switzerland



 QuadRule.w = [2; 2; 2; 4.5; 4.5 ; 4.5;...
         4.5; 4.5; 4.5; 27]/120;
	
 QuadRule.x = [0 0; 1 0; 0 1;...
         1/3 0; 2/3 0; 2/3 1/3;...
         1/3 2/3; 0 2/3; 0 1/3;...
         1/3 1/3];

 return

