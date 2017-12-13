function Mloc = STIMA_Vol_Conv_hpDG_1D(Vertices,p,QuadRule,Shap,GradShap,Vhandle)
% STIMA_Vol_Conv_HPDG_1D
%
%   MLOC = STIMA_Vol_Conv_HPDG_1D(VERTICES,P,QUADRULE,SHAP)

%   Copyright 2007-2007 Patrick Meury & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory

  Mloc = zeros(p+1,p+1);

  % Compute entries of element mass matrix
  
  h = (Vertices(2)-Vertices(1))/2;
  x = (Vertices(2)+Vertices(1))/2 + h*QuadRule.x;
  v=Vhandle(x);
  for i = 1:(p+1)
    for j = 1:(p+1)
      Mloc(i,j) = sum(QuadRule.w.*v.*GradShap(:,i).*Shap(:,j));
    end
  end
  
return