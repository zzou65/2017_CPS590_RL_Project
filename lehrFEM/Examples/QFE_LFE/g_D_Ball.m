function g_D = g_D_Ball(x,varargin)
% G_D_BALL Dirichlet Boundary Data.
%
%   G_D = G_D_BALL(X) computes the value of the Dirichlet boundary data for
%   the function 
%
%     u(r,theta) = cos(pi/2*r)
%
%   on the unit ball.
%
%   Example:
%
%   g_D = g_D_Ball([0 0]);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
    
  r = sqrt(x(:,1).^2+x(:,2).^2);
  g_D = cos(pi/2*r);

return