function grad_shap = grad_shap_QFE(x)
% GRAD_SHAP_QFE Gradient of shape functions.
%
%   GRAD_SHAP = GRAD_SHAP_QFE(X) computes the values of the gradient of the 
%   shape functions for the Lagrangian finite element of order 2 at the
%   quadrature points X.
%
%   Example:
%
%   grad_shap = grad_shap_QFE([0 0]);
%
%   See also shap_QFE.

%   Copyright 2005-2005 Patrick Meury and Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  grad_shap = zeros(size(x,1),12);
  
  grad_shap(:,1) = 4*x(:,1)+4*x(:,2)-3;
  grad_shap(:,2) = 4*x(:,1)+4*x(:,2)-3;
  grad_shap(:,3) = 4*x(:,1)-1;
  grad_shap(:,6) = 4*x(:,2)-1;
  grad_shap(:,7) = 4*(1-2*x(:,1)-x(:,2)); 
  grad_shap(:,8) = (-4)*x(:,1);
  grad_shap(:,9) = 4*x(:,2);
  grad_shap(:,10) = 4*x(:,1);
  grad_shap(:,11) = (-4)*x(:,2);
  grad_shap(:,12) = 4*(1-x(:,1)-2*x(:,2));
  
return
