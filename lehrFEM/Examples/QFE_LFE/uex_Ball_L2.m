function uex = uex_Ball_L2(x,varargin)
% UEX_BALL_L2 Exact solution for L2 norm.
%
%   UEX = UEX_BALL_L2(X) computes the value of the function
%
%      u(r,theta) = cos(pi/2*r)
%
%   at the point X.
%
%   Example:
%
%   uex = uex_Ball_L2([0 0]);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 
  
  % Compute the exact solution for L2 norm
  
  r = sqrt(x(:,1).^2 + x(:,2).^2);
  uex = cos(pi*r/2);

return