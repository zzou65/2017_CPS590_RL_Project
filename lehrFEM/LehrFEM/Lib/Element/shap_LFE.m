function shap = shap_LFE(x)
% SHAP_LFE Shape functions.
%
%   SHAP = SHAP_LFE(X) computes the values of the shape functions for
%   the Lagrangian finite element of order 1 at the quadrature points X.
%
%   Example:
%
%   shap = shap_LFE([0 0]);
%
%   See also grad_shap_LFE.    

%   Copyright 2005-2005 Patrick Meury and Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  shap = zeros(size(x,1),3);

  shap(:,1) = 1-x(:,1)-x(:,2);
  shap(:,2) = x(:,1);
  shap(:,3) = x(:,2);

return