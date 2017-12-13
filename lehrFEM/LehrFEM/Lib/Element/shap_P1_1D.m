function shap = shap_P1_1D(x)
% SHAP_1D Shape functions.
%
%   SHAP = SHAP_P1_1D(X) computes the values of the shape functions for
%   linear finite elements at the quadrature points X.
%
%   Example:
%
%   shap = shap_P1_1D(0);
%
%   See also grad_shap_P1_1D.    

%   Copyright 2005-2005 Patrick Meury and Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  nPts = size(x,1);
  
  % Preallocate memory
  
  shap = zeros(nPts,2);
  
  % Compute function values
  
  shap(:,1) = 1-x;
  shap(:,2) = x;

return
