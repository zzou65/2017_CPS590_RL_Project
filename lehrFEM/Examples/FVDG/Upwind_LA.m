function H = Upwind_LA(x,U_L,U_R,Normal,varargin)
% UPWIND_LA Upwind flux for linear advection.
%
%   H = UPWIND_LA(X,U_L,U_R,NORMAL) computes the value of the upwind
%   numerical flux for the linear advection flux function.
%
%   Example:
%
%   H = Upwind_LA([0.5 0.5],U_R,U_L,Normal);
%
%   See also f_LA.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  % Initialize constants

  nPts = size(x,1);
  
  % Compute value of upwind numerical flux
 
  dpr = sum(Normal); 
  if(dpr < 0);
    H = f_LA(x,U_L)*transpose(Normal);
  else
    H = f_LA(x,U_R)*transpose(Normal);
  end
  
return
