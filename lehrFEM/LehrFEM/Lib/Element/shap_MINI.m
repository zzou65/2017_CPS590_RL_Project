function shap = shap_MINI(x)
% SHAP_MINI Shape functions.
%
%   SHAP = SHAP_MINI(X) computes the values of the shape functions for
%   the MINI element at the quadrature points X.
%
%   Example:
%
%   shap = shap_MINI([0 0]);
%
%   See also grad_shap_MINI.    

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  shap = zeros(size(x,1),4);

  shap(:,1) = 1-x(:,1)-x(:,2);
  shap(:,2) = x(:,1);
  shap(:,3) = x(:,2);
  shap(:,4) = 27*(1-x(:,1)-x(:,2)).*x(:,1).*x(:,2);

return