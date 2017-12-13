function dphi = dpulse1D(x,varargin)
% DPULSE_1D Derivative of a 1D Wave pulse.
%
%   DPHI = DPULSE_1D(X) computes the derivative of a 1D wave pulse.
%
%   Example:
%
%   dphi = dpulse_1D(x);

%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  XL = 0.75; XR=0.25;        % Left and right hand side boundary 

  % Compute function value

  dAmp_1 = zeros(size(x,1),1);
 
  loc = find(x(find( x< XL))>XR);
  
  dAmp_1(loc) =4*pi*cos(2*pi*(x(loc)-0.25)).*sin(2*pi*(x(loc)-0.25));     
   
  dphi = dAmp_1;
  
return