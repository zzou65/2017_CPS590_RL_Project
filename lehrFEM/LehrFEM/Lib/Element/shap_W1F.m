function shap = shap_W1F(x)
% SHAP_W1F Shape functions.
%
%   SHAP = SHAP_W1F(X) computes the values of the shape functions for the
%   edge finite element (Whitney-1-Form) at the quadrature points X.
%
%   Example:
%
%   shap = shap_W1F([0 0]);
%
%   See also shap_LFE2.    

%   Copyright 2005-2005 Patrick Meury and Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  shap = zeros(size(x,1),6);

  shap(:,1) = -x(:,2);
  shap(:,2) = x(:,1);
  shap(:,3) = -x(:,2);
  shap(:,4) = x(:,1)-1;
  shap(:,5) = 1-x(:,2);
  shap(:,6) = x(:,1);

return
