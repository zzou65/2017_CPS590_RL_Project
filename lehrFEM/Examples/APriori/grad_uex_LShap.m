function grad_uex = grad_uex_LShap_H1S(x,varargin)
% GRAD_UEX_LSHAP_H1S Exact solution for H1 semi-norm.
%
%   GRAD_UEX = GRAD_UEX_LSHAP(X) computes the gradient of the function
%
%      u(r,theta) = r^2/3*sin(2/3*theta)
%
%   at the point X.
%
%   Example:
%
%   grad_uex = grad_uex_LShap([0 0]);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  % Compute polar coordinates

  r = sqrt(x(:,1).^2 + x(:,2).^2);
  theta = atan2(x(:,2),x(:,1));
  Loc = theta < 0;
  theta(Loc) = theta(Loc) + 2*pi;
  
  % Compute the exact solution for H1 semi norm
  
  grad_uex = 2/3*(r.^(-1/3)*ones(1,2)).*[-sin(theta/3) cos(theta/3)];
  
return