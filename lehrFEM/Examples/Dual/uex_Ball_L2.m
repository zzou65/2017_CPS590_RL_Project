function uex = uex_Ball_L2(x,varargin)
% UEX_BALL_L2 Exact solution for L2 norm.
%
%   UEX = UEX_BALL_L2(X) computes the value of the function
%
%     uex=sin(pi*x1)*sinh(*x2)
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
  
  uex = sin(pi.*x(:,1)).*sinh(pi.*x(:,2));

return