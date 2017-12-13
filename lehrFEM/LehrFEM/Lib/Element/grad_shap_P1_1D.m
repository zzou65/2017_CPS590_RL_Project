function grad_shap = grad_shap_P1_1D(x)
% GRAD_SHAP_P1_1D Gradient of shape functions.
%
%   GRAD_SHAP = GRAD_SHAP_P1_1D(X) computes the values of the derivatives
%   of the shape functions for linear finite elements.
%
%   Example:
%
%   grad_shap = grad_shap_P1_1D(0);
%
%   See also shap_P1_1D.

%   Copyright 2005-2005 Patrick Meury and Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
  
  nPts = size(x,1);
  
  % Preallocate memory

  grad_shap = zeros(nPts,2);

  % Compute values of derivatives
  
  grad_shap(:,1) = -ones(nPts,1);
  grad_shap(:,2) = ones(nPts,1);
  
return