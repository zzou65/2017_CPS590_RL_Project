function Mloc = MASS_Vol_hpDG_1D(Vertices,p,QuadRule,Shap,varargin)
% MASS_VOL_HPDG_1D
%
%   MLOC = MASS_VOL_HPDG_1D(VERTICES,P,QUADRULE,SHAP)

%   Copyright 2007-2007 Patrick Meury & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory

  Mloc = zeros(p+1,p+1);

  % Compute entries of element mass matrix
  
  h = (Vertices(2)-Vertices(1))/2;
  for i = 1:(p+1)
    for j = i:(p+1)
      Mloc(i,j) = sum(QuadRule.w.*Shap(:,i).*Shap(:,j))*h;
    end
  end
  
  % Fill in lower triangular part
    
  tri = triu(Mloc);
  Mloc = tril(transpose(tri),-1)+tri;
  
return