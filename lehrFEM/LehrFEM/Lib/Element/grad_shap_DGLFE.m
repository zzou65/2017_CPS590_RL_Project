function grad_shap = grad_shap_DGLFE(x)
% GRAD_SHAP_DGLFE Gradient of shape functions.
%
%   GRAD_SHAP = GRAD_SHAP_DGLFE(X) computes the values of gradients of the
%   shape functions for discontinuous linear Lagrangian finite elements at
%   the quadrature points X.
%
%   Example:
%
%   grad_shap = grad_shap_DGLFE([0 0]);
%
%   See also shap_DGLFE.    

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constants

    nPts = size(x,1);

    % Preallocate memory
  
    grad_shap = zeros(nPts,6);

    % Compute values of gradients
  
    grad_shap(:,1) = -1;
    grad_shap(:,2) = -1;
    grad_shap(:,3) =  1;
    grad_shap(:,6) =  1;
    
return