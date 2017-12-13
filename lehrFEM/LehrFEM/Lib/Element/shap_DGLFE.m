function shap = shap_DGLFE(x)
% SHAP_DGLFE Shape functions.
%
%   SHAP = SHAP_DGLFE(X) computes the values of the shape functions for the
%   discontinuous linear Lagrangian finite element at the quadrature points
%   X.
%
%   Example:
%
%   shap = shap_DGLFE([0 0]);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  shap = zeros(size(x,1),3);

  shap(:,1) = 1-x(:,1)-x(:,2);
  shap(:,2) = x(:,1);
  shap(:,3) = x(:,2);

return