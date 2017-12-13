function [uex,grad_uex] = uex_Ball_H1(x,varargin)
% UEX_BALL_H1 Exact solution for H1 norm.
%
%   [UEX,GRAD_UEX] = UEX_BALL_H1(X) computes the value and the gradient of
%   the function
%
%      u(r,theta) = cos(pi/2*r)
%
%   at the point X.
%
%   Example:
%
%   [uex,grad_uex] = uex_Ball_H1([0 0]);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Compute polar coordinates
  
  r = sqrt(x(:,1).^2 + x(:,2).^2);
  theta = atan2(x(:,2),x(:,1));
  Loc = theta < 0;
  theta(Loc) = theta(Loc) + 2*pi;
  
  % Compute the exact solution for H1 norm
  
  uex = cos(pi*r/2);
  grad_uex = [-pi/2*cos(theta).*sin(pi*r/2) -pi/2*sin(theta).*sin(pi*r/2)];

return