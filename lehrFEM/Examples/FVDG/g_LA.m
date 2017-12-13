function g = g_LA(x,ElemFlag,t,nu,varargin)
% G_LA Load data.
%
%   G = G_LA(X,ELEMFLAG,T,NU) computes the value of the right hand side
%   load data.
%
%   The smaller the real valued parameter NU the steeper the gradient of
%   the exact solution at X(:,1) = -X(:,2). 
%
%   Example:
%
%   g = g_LA([0 0],1,1,0.01);
%
%   See also f_LA.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  g = 4*(1-x(:,1).^2).*x(:,1).*(1-x(:,2).^2).^2.* ...
      atan((x(:,1)+x(:,2)+2*(1-t))/nu) + ...
      4*(1-x(:,1).^2).^2.*(1-x(:,2).^2).*x(:,2).* ...
      atan((x(:,1)+x(:,2)+2*(1-t))/nu);

return