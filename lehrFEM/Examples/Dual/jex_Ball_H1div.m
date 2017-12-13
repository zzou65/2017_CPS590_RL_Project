function grad_uex = jex_Ball_H1div(x,varargin)
% UEX_BALL_H1S Exact solution for H1 semi-norm.
%
%   GRAD_UEX = UEX_BALL_H1S(X) computes the gradient of the function
%
%      u(r,theta) = sin(pi.*x1).*sinh(pi.*x2);
%
%   at the point X.
%
%   Example:
%
%   grad_uex = uex_Ball_H1S([0 0]);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 
  
  grad_uex = [(-pi).*sin(pi.*x(:,1)).*cosh(pi.*x(:,2)) (+pi).*cos(pi.*x(:,1)).*sinh(pi.*x(:,2))];
  
return