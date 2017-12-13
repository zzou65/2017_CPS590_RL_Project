function uex = uex_LShap_L2(x,varargin)
% UEX_LSHAP_L2 Exact solution for L2 norm.
%
%   UEX = UEX_LSHAP_L2(X) computes the value of the function
%
%      u(r,theta) = r^2/3*sin(2/3*theta)
%
%   at the point X.
%
%   Example:
%
%   uex = uex_LShap_L2([0 0]);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 
  
  % Compute polar coordinates

  r = sqrt(x(:,1).^2 + x(:,2).^2);
  theta = atan2(x(:,2),x(:,1));
  Loc = theta < 0;
  theta(Loc) = theta(Loc) + 2*pi;
 
  % Compute the exact solution for L2 norm
  
  uex = r.^(2/3).*sin(2*theta/3);
  
return