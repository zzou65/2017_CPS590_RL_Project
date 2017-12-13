function g = g(x,BdFlag,t,varargin)
% G Dirichlet boundary data.
%
%   G = G(X,BDFLAG,T) computes the value of the Dirichlet boundary data at
%   the points X for the the time T.
%
%   Example:
%
%   g = g(x,-1,0);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  g = zeros(size(x,1),1);
  
return