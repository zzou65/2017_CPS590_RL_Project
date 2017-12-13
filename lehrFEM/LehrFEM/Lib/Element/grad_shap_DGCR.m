function grad_shap = grad_shap_DG(x)
% GRAD_SHAP_DGCR Gradient of shape functions.
%
%   GRAD_SHAP = GRAD_SHAP_DGCR(X) computes the values of gradients of the
%   shape functions for discontinuous Crouzeix-Raviart finite elements at
%   the quadrature points X.
%
%   Example:
%
%   grad_shap = grad_shap_DGCR([0 0]);
%
%   See also shap_DGCR.    

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constants

    nPts = size(x,1);

    % Preallocate memory
  
    grad_shap = zeros(nPts,6);

    % Compute values of gradients
  
    grad_shap(:,1) =  2;
    grad_shap(:,2) =  2;
    grad_shap(:,3) = -2;
    grad_shap(:,6) = -2;
    
return