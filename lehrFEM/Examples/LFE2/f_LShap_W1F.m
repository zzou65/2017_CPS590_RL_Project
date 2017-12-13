function f = f_LShap_W1F(x,varargin)
% F_LSHAP_W1F Exact solution for H1 semi-norm.
%
%   F = F_LSHAP_W1F(X) computes the right hand side in the W1F problem of
%   the function
%
%      u(r,theta) = r^2/3*sin(2/3*theta)
%
%   at the point X.
%
%   Example:
%
%   f = f_LShap_W1F([0 0]);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  % Compute polar coordinates

  r = sqrt(x(:,1).^2 + x(:,2).^2);
  theta = atan2(x(:,2),x(:,1));
  Loc = theta < 0;
  theta(Loc) = theta(Loc) + 2*pi;
  
  % Compute the exact solution for H1 semi norm
  
  f = 2/3*(r.^(-1/3)*ones(1,2)).*[-sin(theta/3) cos(theta/3)];
  
return