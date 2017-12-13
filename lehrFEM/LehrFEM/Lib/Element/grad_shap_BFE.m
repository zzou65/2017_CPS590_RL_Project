function grad_shap = grad_shap_BFE(x)
% GRAD_SHAP_BFE Gradient of shape functions.
%
%   GRAD_SHAP = GRAD_SHAP_BFE(X) computes the values of the gradient of the
%   shape functions for the Bilinear finite elements at the quadrature
%   points X.
%
%   Example:
%
%   grad_shap = grad_shap_BFE([0 0]);
%
%   See also shap_BFE.

%   Copyright 2005-2005 Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(x,1);

  % Preallocate memory
  
  grad_shap = zeros(nPts,8);
  
  % Compute values of gradients
  
  grad_shap(:,1) = x(:,2)-1; 
  grad_shap(:,2) = x(:,1)-1;
  grad_shap(:,3) = 1-x(:,2);
  grad_shap(:,4) = -x(:,1);
  grad_shap(:,5) = x(:,2);
  grad_shap(:,6) = x(:,1);
  grad_shap(:,7) = -x(:,2);
  grad_shap(:,8) = 1-x(:,1);
  
return
