function shap = shap_CR(x)
% SHAP_CR Shape functions.
%
%   SHAP = SHAP_CR(X) computes the values of the shape functions for the
%   Crouzeix-Raviart finite element at the quadrature points X.
%
%   Example:
%
%   shap = shap_CR([0 0]);
%
%   See also grad_shap_CR.    

%   Copyright 2005-2005 Patrick Meury and Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  shap = zeros(size(x,1),3);

  shap(:,1) = 2*(x(:,1)+x(:,2))-1;
  shap(:,2) = 1-2*x(:,1);
  shap(:,3) = 1-2*x(:,2);

return