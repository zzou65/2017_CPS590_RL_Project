function shap = shap_LFE2(x)
% SHAP_LFE2 Shape functions.
%
%   SHAP = SHAP_LFE2(X) computes the values of the shape functions for
%   the vector valued Lagrangian finite element of order 1 at the 
%   quadrature points X.
%
%   Example:
%
%   shap = shap_LFE2([0 0]);
%
%   See also shap_LFE, shap_W1F.    

%   Copyright 2005-2005 Patrick Meury and Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  shap = zeros(size(x,1),12);

  shap(:,1) = 1-x(:,1)-x(:,2);
  shap(:,4) = 1-x(:,1)-x(:,2);
  shap(:,5) = x(:,1);
  shap(:,8) = x(:,1);
  shap(:,9) = x(:,2);
  shap(:,12) = x(:,2);

return
