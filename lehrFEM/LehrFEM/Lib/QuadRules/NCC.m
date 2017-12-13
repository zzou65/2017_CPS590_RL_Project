function QuadRule = NCC(order)
% wrap for ncc_trangle_rule of Burkhard
%
%   several Newton-Coates-type  quadrature rules
%   for polynomial with degree order-1.
% 
%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  [x w] = ncc_triangle_rule(order,10);
  QuadRule.x=x';
  QuadRule.w=w'/2;
            
return