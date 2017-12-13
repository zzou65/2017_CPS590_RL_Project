function grad_shap = grad_shap_EP2(x)
% GRAD_SHAP_EP2 Gradient of shape functions.
%
%   GRAD_SHAP = GRAD_SHAP_EP2(X) computes the values of the gradient of the 
%   shape functions for the Lagrangian finite element of order 2 at the
%   quadrature points X.
%
%   Example:
%
%   grad_shap = grad_shap_EP2([0 0]);
%
%   See also shap_EP2.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  grad_shap = zeros(size(x,1),6);
  
  grad_shap(:,1) = 4*x(:,2); 
  grad_shap(:,2) = 4*x(:,1);
  grad_shap(:,3) = (-4)*x(:,2);
  grad_shap(:,4) = 4*(1-x(:,1)-2*x(:,2));
  grad_shap(:,5) = 4*(1-2*x(:,1)-x(:,2));
  grad_shap(:,6) = (-4)*x(:,1);
  
return