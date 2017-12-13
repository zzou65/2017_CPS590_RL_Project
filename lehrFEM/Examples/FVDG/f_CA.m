function f = f_CA(x,u,varargin) 
% F_CA Circular advection flux.
%
%   F = F_CA(X,U) computes the value of the 2D circular advection flux.
%
%   Example:
%
%   f = f_CA([0.5 0.5],0);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  % Initialize constants
  
  OMEGA = 2*pi;  % Angular velocity

  % Compute value of velocity field
  
  nPts = size(x,1);
  amp = ones(nPts,1);
  vfield = zeros(nPts,2);
  vfield(:,1) = -x(:,2);
  vfield(:,2) =  x(:,1);
  
  % Compute value of circular flux
  
  f = zeros(nPts,2);
  f(:,1) = OMEGA*u.*vfield(:,1);
  f(:,2) = OMEGA*u.*vfield(:,2);
  
return