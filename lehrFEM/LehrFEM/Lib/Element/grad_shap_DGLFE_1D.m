function grad_shap = grad_shap_DGLFE_1D(x)
% GRAD_SHAP_DGLFE Gradient of shape functions.
%
%   GRAD_SHAP = GRAD_SHAP_DGLFE(X) computes the values of gradients of the
%   shape functions for discontinuous linear Lagrangian finite elements at
%   the quadrature points X. [-1,1] as reference element.
%
%   Example:
%
%   grad_shap = grad_shap_DGLFE_1D([0;1]);
%
%   See also shap_DGLFE.    

%   Copyright 2007-2007 Patrick Meury Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constants

    nPts = size(x,1);

    % Preallocate memory
  
    grad_shap = zeros(nPts,2);

    % Compute values of gradients
  
    grad_shap(:,1) = -1/2;
    grad_shap(:,2) = 1/2;
        
return