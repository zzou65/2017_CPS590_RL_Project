function shap = shap_DGLFE_1D(x,varargin)
% SHAP_DGLFE_1D Shape functions.
%
%   SHAP = SHAP_DGLFE_1D(X) computes the values of the shape functions for the
%   discontinuous linear Lagrangian finite element at the quadrature points
%   X. [-1,1] as reference element.
%
%   Example:
%
%   shap = shap_DGLFE_1D([0; 1]);

%   Copyright 2007-2007 Patrick Meury Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  shap = zeros(size(x,1),2);

  shap(:,1) = 1/2*(-x+1);
  shap(:,2) = 1/2*(x+1);

return