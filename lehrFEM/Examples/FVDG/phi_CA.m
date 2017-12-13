function phi = phi_CA(z,r,varargin)
% PHI_CA Compactly supported blob-like shaped function.
%
%   PHI = PHI_CA(Z,R) computes the value of a blob-like shaped compactly
%   supported function.
%
%   Example:
%
%   phi = phi_CA(z,0.4);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  phi = zeros(size(z));
  Loc = find(z <= r);
  
  %phi(Loc) = cos(pi/2*z(Loc)/r).^2;  % Smooth blob shaped profile
  phi(Loc) = 1;                    % Non-smooth cylindrical profile
  
return

