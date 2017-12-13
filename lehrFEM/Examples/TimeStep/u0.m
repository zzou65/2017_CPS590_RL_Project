function u0 = u0(x,ElemFlag,varargin)
% U0 Initial data.
%
%   U0 = U0(X,ELEMFLAG) computes the value of the initial data at the
%   points X.
%
%   Example:
%
%   u0 = u0(x,-1);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  u0 = sin(2*pi*x(:,1)).*sin(2*pi*x(:,2));
  
return