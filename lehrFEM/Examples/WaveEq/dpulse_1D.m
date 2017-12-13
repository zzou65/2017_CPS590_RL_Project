function dphi = dpulse_1D(x,varargin)
% DPULSE_1D Derivative of a 1D Wave pulse.
%
%   DPHI = DPULSE_1D(X) computes the derivative of a 1D wave pulse.
%
%   Example:
%
%   dphi = dpulse_1D(x);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  X0 = 1  ;      % Left and right hand side boundary 
  S = 1;         % Scaling parameter
  OMEGA = 8*pi;  % Frequency
  
  % Compute function value

  Amp_1 = zeros(size(x,1),1);
  dAmp_1 = zeros(size(x,1),1);
  Amp_2 = zeros(size(x,1),1);
  dAmp_2 = zeros(size(x,1),1);
  
  loc = find(abs(x) < X0);
  
  Amp_1(loc) = ((X0^2-x(loc).^2)/X0^2).^3;
  dAmp_1(loc) = -6/X0^2*x(loc).*((X0^2-x(loc).^2)/X0^2).^2;
  
  Amp_2(loc) = exp(-S*x(loc).^2);
  dAmp_2(loc) = -2*S*x(loc).*exp(-S*x(loc).^2);
 
  dphi = dAmp_1.*Amp_2 + Amp_1.*dAmp_2;
  
return