 function Mloc = MASS_QFELFE(Vertices,n,QuadRule)
% MASS_QFELFE Element mass matrix.
%
%   MLOC = MASS_QFELFE(VERTICES,n,Quadrule) computes the element mass matrix using
%   quadratic and linear Lagrangian finite elements .
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current
%   element in a row wise orientation.
%
%   Example:
%
%   Mloc = MASS_QFELFE(Vertices,1,QuadRule);

%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

% local MASS-Matrix using Quadrature rule

  Mloc=zeros(3,6);

  s=size(QuadRule.x,1);
  n=shap_LFE(QuadRule.x);
  for j=1:s
         Mloc=Mloc+QuadRule.w(j)*n(j,:)'*[n(j,1)*(2*n(j,1)-1) ...
             n(j,2)*(2*n(j,2)-1) ...
             n(j,3)*(2*n(j,3)-1) ...
             4*n(j,1)*n(j,2) 4*n(j,2)*n(j,3) 4*n(j,1)*n(j,3)];
  end