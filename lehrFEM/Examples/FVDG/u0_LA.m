function u = u0_LA(x,ElemFlag,nu,varargin)
% U0_LA Intial data.
%
%   U = U0_LA(X,ELEMFLAG,NU) computes the value of the initial data.
%
%   The smaller the real valued parameter NU the steeper the gradient of
%   the exact solution at X(:,1) = -X(:,2). 
%
%   Example:
%
%   u = u0_LA([0 0],0,0.01);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  u = -(1-x(:,1).^2).^2.*(1-x(:,2).^2).^2.*atan((x(:,1)+x(:,2)+2)/nu);

return