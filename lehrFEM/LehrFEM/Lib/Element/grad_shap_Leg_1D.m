function grad_shap = grad_shap_Leg_1D(x,p)
% GRAD_SHAP_LEG_1D Gradients of shape functions for 1D hpDG discretization.
%
%   GRAD_SHAP = GRAD_SHAP_LEG_1D(X,P) gradients of shape functions for 1D
%   hpDG discretizations using Legedre polynomials up to degree P as basis
%   functions.
%
%   Example:
%
%   grad_shap = grad_shap_Leg_1D(0,10);
%
%   See also shap_Leg_1D.

%   Copyright 2007-2007 Patrick Meury & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(x,1);
  
  % Compute Legendre polynomials
  
  P0 = ones(nPts,1);
  P1 = x;

  % Preallocate memory
  
  grad_shap = zeros(nPts,p+1);
  
  % Compute gradients of shape functions
  
  grad_shap(:,1) = zeros(nPts,1);
  
  if (p>0)
   grad_shap(:,2) = ones(nPts,1);
  
   for i = 1:(p-1)
    grad_shap(:,i+2) = (2*i+1)*P1 + grad_shap(:,i);
    tmp = P1;
    P1 = ((2*i+1)*x.*P1 - i*P0)/(i+1);
    P0 = tmp;
   end
  end
   
return