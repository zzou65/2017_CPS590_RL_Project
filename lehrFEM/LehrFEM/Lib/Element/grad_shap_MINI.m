function grad_shap = grad_shap_MINI(x)
% GRAD_SHAP_MINI Gradient of shape functions.
%
%   GRAD_SHAP = GRAD_SHAP_MINI(X) computes the values of the gradient
%   of the shape functions for the MINI element at the quadrature points X.
%
%   Example:
%
%   grad_shap = grad_shap_MINI([0 0]);
%
%   See also shap_MINI.

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(x,1);

  % Preallocate memory
  
  grad_shap = zeros(nPts,8);
  
  % Compute values of gradients
  
  grad_shap(:,1:2) = -ones(nPts,2);
  grad_shap(:,3) = ones(nPts,1);
  grad_shap(:,6) = ones(nPts,1);
  grad_shap(:,7) = 27*(x(:,2)-2*x(:,1).*x(:,2)-x(:,2).^2);
  grad_shap(:,8) = 27*(x(:,1)-x(:,1).^2-2*x(:,1).*x(:,2));

return
