function shap = shap_DGCR(x)
% SHAP_DGCR Shape functions.
%
%   SHAP = SHAP_DGCR(X) computes the values of the shape functions for the
%   discontinuous Crouzeix-Raviart finite element at the quadrature points
%   X.
%
%   Example:
%
%   shap = shap_DGCR([0 0]);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  shap = zeros(size(x,1),3);

  shap(:,1) = 2*(x(:,1)+x(:,2))-1;
  shap(:,2) = 1-2*x(:,1);
  shap(:,3) = 1-2*x(:,2);

return