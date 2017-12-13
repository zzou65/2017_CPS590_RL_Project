function QuadRule = TProd(QuadRule)
% TPROD 2D Tensor-product Gauss quadrature.
%
%   QUADRULE = TPROD(QUADRULE) computes a 2D tensor-product Gauss quadrature
%   strating from the 1D Gauss quadrature QUADRULE.
%
%   The struct QUADRULE should at least contain the following fields:
%    W N-by-1 matrix specifying the weights of a 1D quadrature rule.
%    X N-by-1 matrix specifying the abscissae of a 1D quadrature rule.
%
%   Example:
%
%   QuadRule = TProd(gauleg(10));

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  nGauss_1D = size(QuadRule.w,1);
  nGauss_2D = nGauss_1D^2;
    
  % Preallocate memory
  
  w = zeros(nGauss_2D,1);
  x = zeros(nGauss_2D,2);
  
  % Build 2D quadrature on the reference square
  
  k = 1;
  for i = 1:nGauss_1D
    for j = 1:nGauss_1D
      w(k) = QuadRule.w(i)*QuadRule.w(j);
      x(k,1) = QuadRule.x(i);
      x(k,2) = QuadRule.x(j);
      k = k+1;
    end
  end
  
  % Assign output arguements
  
  QuadRule.w = w;
  QuadRule.x = x;
  
return