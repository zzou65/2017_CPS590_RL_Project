function Mloc = MASS_Vol_PDG(Vertices,ElemFlag,QuadRule,Shap,varargin)
% MASS_VOL_PDG Element mass matrix.
%
%   MLOC = MASS_VOL_PDG(VERTICES,ELEMFLAG,QUADRULE,SHAP) computes the
%   element mass matrix using the Quadrature rule QUADRULE and the shape
%   functions SHAP.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current
%   element in a row wise orientation.
%
%   The integer ELEMFLAG denotes the element flag of the current element.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%
%   SHAP is a function handle to the reference element shape functions.
%
%   Example:
%   
%   Mloc = MASS_Vol_PDG(Vertices,ElemFlag,P3O3(),@shap_DGCR);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Compute element mapping

  bK = Vertices(1,:);
  BK = [Vertices(2,:) - bK; ...
        Vertices(3,:) - bK];
  det_BK = abs(det(BK));
  
  % Evaluate shape functions
  
  N = Shap(QuadRule.x);
  nDofs = size(N,2);
  
  % Preallocate memory
  
  Mloc = zeros(nDofs,nDofs);
  
  % Compute local mass matrix
  
  for j1 = 1:nDofs
    for j2 = j1:nDofs
      Mloc(j1,j2) = sum(QuadRule.w.*N(:,j1).*N(:,j2))*det_BK;
    end
  end

  % Fill in lower trinagular part
    
  tri = triu(Mloc);
  Mloc = tril(transpose(tri),-1)+tri;
  
return