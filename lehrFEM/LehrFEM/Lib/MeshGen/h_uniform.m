function h = h_uniform(x,varargin)
% H_UNIFORM Element size function.
%
%   H = H_UNIFORM(X,VARARGIN) computes the value of the element size
%   function for uniform relative distribution.
%
%   Example:
%
%   h = h_unifrom(x);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  h = ones(size(x,1),1);
  
return