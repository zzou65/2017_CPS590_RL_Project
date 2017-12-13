function grad_shap = grad_shap_LFE(x)
% GRAD_SHAP_LFE Gradient of shape functions.
%
%   GRAD_SHAP = GRAD_SHAP_LFE(X) computes the values of the gradient
%   of the shape functions for the Lagrangian finite element of order 1
%   at the quadrature points X.
%
%   Example:
%
%   grad_shap = grad_shap_LFE([0 0]);
%
%   See also shap_LFE.

%   Copyright 2005-2005 Patrick Meury and Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(x,1);

  % Preallocate memory
  
  grad_shap = zeros(nPts,6);
  
  % Compute values of gradients
  
  grad_shap(:,1:2) = -ones(nPts,2);
  grad_shap(:,3) = ones(nPts,1);
  grad_shap(:,6) = ones(nPts,1);
  
return
