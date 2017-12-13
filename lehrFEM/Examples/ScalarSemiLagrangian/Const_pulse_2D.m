function phi = Const_pulse_2D(x,R0,C,varargin)
% PULSE_2D Value of a 2D pulse.
%
%   PHI = PULSE_2D(X,R0,C) computes the value of a 2D wave pulse.
%
%   Example:
%
%   phi = pulse_2D(x);

%   Copyright 2006-2010 Patrick Meury, Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory
  
  phi = zeros(size(x,1),1);
    
  r = sqrt((x(:,1)-C(1)).^2+(x(:,2)-C(2)).^2);
  Loc = find(r < R0);
  phi(Loc) = 1;
  
return