function shap = shap_EP2(x)
% SHAP_EP2 Shape functions.
%
%   SHAP = SHAP_EP2(X) computes the values of the edge shape functions for
%   the Lagrangian finite element of order 2 at the quadrature points X.
%
%   Example:
%
%   shap = shap_EP2([0 0]);
%
%   See also grad_shap_EP2.    

%   Copyright 2005-2005 Patrick Meury 
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  shap = zeros(size(x,1),3);

  shap(:,1) = 4*x(:,1).*x(:,2);
  shap(:,2) = 4*x(:,2).*(1-x(:,1)-x(:,2));
  shap(:,3) = 4*x(:,1).*(1-x(:,1)-x(:,2));

return