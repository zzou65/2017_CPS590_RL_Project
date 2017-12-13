function phi = pulse_1D(x,varargin)
% PULSE_1D Value of a 1D Wave pulse.
%
%   PHI = PULSE_1D(X) computes the value of a 1D wave pulse.
%
%   Example:
%
%   phi = pulse_1D(x);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  X0 = 1;        % Left and right hand side boundary 
  S = 1;         % Scaling parameter
  OMEGA = 8*pi;  % Frequency
  
  % Compute function value

  Amp_1 = zeros(size(x,1),1);
  Amp_2 = zeros(size(x,1),1);
 
  loc = find(abs(x) < X0);
  
  Amp_1(loc) = ((X0^2-x(loc).^2)/X0^2).^3;
  Amp_2(loc) = exp(-S*x(loc).^2);
     
  phi = Amp_1.*Amp_2;
  
return