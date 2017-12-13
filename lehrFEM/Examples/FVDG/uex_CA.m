function u = u0_CA(x,ElemFlag,t,varargin)
% UEX_CA Exact solution.
%
%   U = UEX_CA(X,ELEMFLAG,T) computes the value of the exact solution.
%
%   Example:
%
%   u = uex_CA([0 0],0,1.123);
%
%   See also phi_CA.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  % Initialize constants
 
  R0 = 0.5;      % Initial radius
  PHI = pi/4;    % Initial center
  R1 = 0.4;      % Blob radius (support)
  OMEGA = 2*pi;  % Angular velocity
  
  % Compute relative coordinates
  
  C = R0*[cos(OMEGA*t+PHI) sin(OMEGA*t+PHI)];
  r = sqrt((x(:,1)-C(1)).^2+(x(:,2)-C(2)).^2);
  
  % Compute value of initial blob
  
  u = phi_CA(r,R1);
  
return