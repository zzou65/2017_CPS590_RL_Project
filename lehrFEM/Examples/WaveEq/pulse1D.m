function phi = pulse1D(x,varargin)
% PULSE_1D Value of a 1D Wave pulse.
%
%   PHI = PULSE_1D(X) computes the value of a 1D wave pulse.
%
%   Example:
%
%   phi = pulse_1D(x);

%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  XL = 0.75; XR=0.25;        % Left and right hand side boundary 

  % Compute function value

  Amp_1 = zeros(size(x,1),1);
 
  loc = find(x(find( x< XL))>XR);
  
  Amp_1(loc) =1-(cos(2*pi*(x(loc)-0.25))).^2;     
  
  phi = Amp_1;
  
return