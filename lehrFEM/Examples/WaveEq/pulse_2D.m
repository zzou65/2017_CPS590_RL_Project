function phi = pulse_2D(x,varargin)
% PULSE_2D Value of a 2D Wave pulse.
%
%   PHI = PULSE_2D(X) computes the value of a 2D wave pulse.
%
%   Example:
%
%   phi = pulse_2D(x);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  R0 = 0.5;       % Radius
  A = 10;         % Amplitude 
  C = [0.2 0.2];  % Center
  S = 0.0001;     % Scaling parameter for Gaussian
  
  % Preallocate memory
  
  phi = zeros(size(x,1),1);
    
  r = sqrt((x(:,1)-C(1)).^2+(x(:,2)-C(2)).^2);
  Loc = find(r < R0);
  phi(Loc) = A*(1-r(Loc)/R0).^3.*exp(-S*r(Loc).^2);
  
return