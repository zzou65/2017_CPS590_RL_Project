function grad_shap = grad_shap_CR(x)
% GRAD_SHAP_CR Gradient of shape functions.
%
%   GRAD_SHAP = GRAD_SHAP_CR(X) computes the values of gradients of the
%   shape functions for Crouzeix-Raviart finite elements at the quadrature
%   points X.
%
%   Example:
%
%   grad_shap = grad_shap_CR([0 0]);
%
%   See also shap_CR.    

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  nPts = size(x,1);

  % Preallocate memory
  
  grad_shap = zeros(nPts,6);

  % Compute values of gradients
  
  grad_shap(:,1:2) = 2*ones(nPts,2);
  grad_shap(:,3) = -2*ones(nPts,1);
  grad_shap(:,6) = -2*ones(nPts,1);

return