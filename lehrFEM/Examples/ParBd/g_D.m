function g_D = g_D(x,varargin)
% G_D Dirichlet boundary data.
%
%   G_D = G_D(X) computes the value of the Dirichlet boundary data for
%   the function 
%
%     u(r,theta) = cos(pi/2*r)
%
%   on the unit ball.
%
%   Example:
%
%   g_D = g_D([0 0]);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
    
  g_D = zeros(size(x,1),1);

return