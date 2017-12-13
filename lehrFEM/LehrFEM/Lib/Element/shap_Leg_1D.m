function shap = shap_Leg_1D(x,p)
% SHAP_LEG_1D Shape functions for 1D hpDG discretization.
%
%   SHAP = SHAP_LEG_1D(X,P) shape functions for 1D hpDG discretizations
%   using Legedre polynomials up to degree P as basis functions. [-1,1] is
%   reference element.
%
%   Example:
%
%   shap = shap_Leg_1D(0,10);
%
%   See also grad_shap_Leg_1D.

%   Copyright 2007-2007 Patrick Meury & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(x,1);
  
  % Preallocate memory

  shap = zeros(nPts,p+1);

  % Compute shape functions
  
  shap(:,1) = ones(nPts,1);
  if(p > 0)
    shap(:,2) = x;
    for i = 2:1:p
      shap(:,i+1) = ((2*i-1)*x.*shap(:,i) - (i-1)*shap(:,i-1))/i;  
    end
  end
  
return
