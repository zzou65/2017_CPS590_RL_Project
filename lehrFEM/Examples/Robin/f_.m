function f = f_(x,varargin)
% F_ Right hand side source.
%
%   F = F_(X) returns the value of the right hand side source term of
%   the Laplace equation.  
%
%   Example:
%
%   f = f_([0 0]);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
% i)   
%f=ones(size(x,1),1);
% ii)
r = sqrt(x(:,1).^2 + x(:,2).^2);
f = (2-r.^2)./(1+r.^(2)).^(5/2);
return