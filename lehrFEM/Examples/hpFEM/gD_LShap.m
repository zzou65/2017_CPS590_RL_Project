function gD = gD_LShap(x,varargin)
% GD_LSHAP Dirichlet Boundary data.
%
%   GD = GD_LSHAP(X) computes the value of the Dirichlet boundary data for
%   the function 
%
%     u(r,theta) = r^2/3*sin(2/3*theta)
%
%   on the L-shaped domain.
%
%   Example:
%
%   gD = gD_LShap([0 0]);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute polar coordinates

  r = sqrt(x(:,1).^2 + x(:,2).^2);
  theta = atan2(x(:,2),x(:,1));
  Loc = theta < 0;
  theta(Loc) = theta(Loc) + 2*pi;
  
  % Compute Dirichlet data
  
  gD = (r.^(2/3)).*sin(2/3*theta);

return