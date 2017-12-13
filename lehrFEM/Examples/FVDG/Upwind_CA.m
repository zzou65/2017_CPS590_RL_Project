function H = Upwind_LA(x,U_L,U_R,Normal,varargin)
% UPWIND_CA Upwind flux for circular advection.
%
%   H = UPWIND_CA(X,U_L,U_R,NORMAL) computes the value of the upwind
%   numerical flux for the circular advection flux function.
%
%   Example:
%
%   H = Upwind_LA([0.5 0.5],U_R,U_L,Normal);
%
%   See also f_CA.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 
  
  % Compute value of upwind numerical flux
  
  H = zeros(size(x,1),1);
  
  dpr = x(:,1)*Normal(2)-x(:,2)*Normal(1);
  
  Loc = find(dpr < 0);
  H(Loc) = f_CA(x(Loc,:),U_L(Loc))*transpose(Normal);
  
  Loc = find(dpr > 0);
  H(Loc) = f_CA(x(Loc,:),U_R(Loc))*transpose(Normal);
  
return
