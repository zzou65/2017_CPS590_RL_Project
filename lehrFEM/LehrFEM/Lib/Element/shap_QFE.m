function shap = shap_QFE(x)
% SHAP_QFE Shape functions.
%
%   SHAP = SHAP_QFE(X) computes the values of the shape functions for the
%   Lagrangian finite element of order 2 at the quadrature points X.
%
%   Example:
%
%   shap = shap_QFE([0 0]);
%
%   See also grad_shap_QFE.    

%   Copyright 2005-2005 Patrick Meury and Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  shap = zeros(size(x,1),6);

  shap(:,1) = (1-x(:,1)-x(:,2)).*(1-2*x(:,1)-2*x(:,2));
  shap(:,2) = x(:,1).*(2*x(:,1)-1);
  shap(:,3) = x(:,2).*(2*x(:,2)-1);
  shap(:,4) = 4*x(:,1).*(1-x(:,1)-x(:,2));
  shap(:,5) = 4*x(:,1).*x(:,2);
  shap(:,6) = 4*x(:,2).*(1-x(:,1)-x(:,2));

return