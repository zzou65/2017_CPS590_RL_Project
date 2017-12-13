function f = f_LA(x,u,varargin)
% F_LA Linear advection flux.
%
%   F = F_LA(X,U) computes the value of the 2D linear advection flux.
%
%   Example:
%
%   f = f_LA([0.5 0.5],0);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  f = u*ones(1,2);

return