function shap = shap_BFE(x)
% SHAP_BFE Shape functions.
%
%   SHAP = SHAP_BFE(X) computes the values of the shape functions for the
%   bilinear finite element at the quadrature points X.
%
%   Example:
%
%   shap = shap_BFE([0 0]);
%
%   See also grad_shap_BFE.    

%   Copyright 2005-2005 Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  shap = zeros(size(x,1),4);

  shap(:,1) = (1-x(:,1)).*(1-x(:,2));
  shap(:,2) = x(:,1).*(1-x(:,2));
  shap(:,3) = x(:,1).*x(:,2);
  shap(:,4) = (1-x(:,1)).*x(:,2);

return