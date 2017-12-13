function uex = jex_LShap_H1div(x,varargin)
% UEX_LSHAP_H1S Exact solution for H1 semi-norm.
%
%   UEX = JEX_LSHAP_H1DIV(X) computes the gradient of the function
%
%      u(r,theta) = r^2/3*sin(2/3*theta)
%
%   at the point X.
%
%   Example:
%
%   uex = jex_LShap_H1div([0 0]);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  % Compute polar coordinates

  r = sqrt(x(:,1).^2 + x(:,2).^2);
  theta = atan2(x(:,2),x(:,1));
  Loc = theta < 0;
  theta(Loc) = theta(Loc) + 2*pi;
  
  uex = 2/3*(r.^(-1/3)*ones(1,2)).*[-cos(theta/3), -sin(theta/3)];
  
return