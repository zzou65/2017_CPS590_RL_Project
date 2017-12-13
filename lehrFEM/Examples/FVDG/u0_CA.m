function u = u0_CA(x,varargin)
% U0_CA Intial data.
%
%   U = U0_CA(X) computes the value of the initial data.
%
%   Example:
%
%   u = u0_CA([0 0]);
%
%   See also phi_CA.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  % Initialize constants
  
  R0 = 0.5;    % Initial radius
  R1 = 0.4;    % Blob radius (support)
  PHI = pi/4;  % Initial angle
  
  % Compute relative coordinates
  
  C = R0*[cos(PHI) sin(PHI)];
  r = sqrt((x(:,1)-C(1)).^2+(x(:,2)-C(2)).^2);
  
  % Compute value of initial blob
  
  u = phi_CA(r,R1);
  
return