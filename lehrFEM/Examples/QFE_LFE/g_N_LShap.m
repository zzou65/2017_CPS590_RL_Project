function g_N = g_N_LShap(x,BdFlag,varargin)
% G_N_LSHAP Neumann Boundary Data.
%
%   G_N = G_N_LSHAP(X) computes the value of the Neumann boundary data for
%   the function 
%
%     u(r,theta) = r^2/3*sin(2/3*theta)
%
%   on the L-shaped domain according to the boundary flag BDFLAG.
%
%   Example:
%
%   g_N = g_N_LShap([0 0],-1);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute polar coordinates

  r = sqrt(x(:,1).^2 + x(:,2).^2);
  theta = atan2(x(:,2),x(:,1));
  Loc = theta < 0;
  theta(Loc) = theta(Loc) + 2*pi;
  
  % Compute Neumann data
  
  switch(BdFlag)
    case -1
      g_N = -(2/3)*r.^(-1/3).*cos(theta/3);
    case -2
      g_N = (2/3)*r.^(-1/3).*sin(theta/3);
    case -3
      g_N = -(2/3)*r.^(-1/3).*sin(theta/3);
    case -4
      g_N = (2/3)*r.^(-1/3).*cos(theta/3);
    otherwise
      g_N = zeros(size(r));  
  end
 
return