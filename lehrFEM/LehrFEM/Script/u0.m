function u = u0(x,varargin)
% U0(X) initial value for discontinuous problem
%
%   U = U0(X) computes the discontinuous function
%
%   at the point X.
%
%   Example:
%
%   u = u0([0 0]);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    u = zeros(length(x),1);
    
    loc = find(x(:,1)>.25 & x(:,1)<.75 & x(:,2)>.25 & x(:,2)<.75);
    u(loc) = 1;
    
return