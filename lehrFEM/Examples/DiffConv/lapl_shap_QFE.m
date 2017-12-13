function lapl_shap = lapl_shap_QFE(x)
% LAPL_SHAP_QFE Gradient of shape functions.
%
%   LAPL_SHAP = LAPL_SHAP_QFE(X) computes the values of the Laplacian of the 
%   shape functions for the Lagrangian finite element of order 2 at the
%   quadrature points X.
%
%   Example:
%
%   grad_shap = grad_shap_QFE([0 0]);
%
%   See also shap_QFE.

%   Copyright 2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  grad_shap = zeros(size(x,1),24);
  nPts = size(x,1);
  lapl_shap(:,1) = 4*ones(nPts,1);
  lapl_shap(:,2) = 4*ones(nPts,1);
  lapl_shap(:,3) = 4*ones(nPts,1);
  lapl_shap(:,4) = 4*ones(nPts,1);
  
  lapl_shap(:,5) = 4*ones(nPts,1);
  %lapl_shap(:,6) = 0;
  %lapl_shap(:,7) = 0;
  %lapl_shap(:,8) = 0;
  
  %lapl_shap(:,9) = 0;
  %lapl_shap(:,10) = 0;
  %lapl_shap(:,11) = 0;
  lapl_shap(:,12) = 4*ones(nPts,1);
  
  lapl_shap(:,13) = -8*ones(nPts,1);
  lapl_shap(:,14) = -4*ones(nPts,1);
  lapl_shap(:,15) = -4*ones(nPts,1);
  %lapl_shap(:,16) = 0;
  
  %lapl_shap(:,17) = 0;
  lapl_shap(:,18) = 4*ones(nPts,1);
  lapl_shap(:,19) = 4*ones(nPts,1);
  %lapl_shap(:,20) = 0;
  
  %lapl_shap(:,21) = 0;
  lapl_shap(:,22) = -4*ones(nPts,1);
  lapl_shap(:,23) = -4*ones(nPts,1);
  lapl_shap(:,24) = -8*ones(nPts,1);
  
return