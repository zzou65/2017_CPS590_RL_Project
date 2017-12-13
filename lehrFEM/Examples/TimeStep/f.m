function f = f(x,ElemFlag,t,varargin)
% F Right hand side source term.
%
%   F = F(X,ELEMFLAG,T) computes the value of the right hand side source
%   term at the points X for the the time T.
%
%   Example:
%
%   f = f(x,-1,0);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  f = -2*pi*sin(2*pi*t).*sin(2*pi*x(:,1)).*sin(2*pi*x(:,2)) + ...  
       8*pi^2*cos(2*pi*t).*sin(2*pi*x(:,1)).*sin(2*pi*x(:,2));

return