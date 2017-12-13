function sigma = sigma_y_const(x,ElemFlag,sigma_0,L,varargin)
% SIGMA_Y_CONST Absorption profile.
%
%   SIGMA = SIGMA_Y_CONST(X,ELEMFLAG,SIGMA_0,L) computes the value of the
%   constant absorption profile in x-direction.
%
%   Example:
%
%   sigma = sigma_y_const(x,0,sigma_0,L);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  sigma = zeros(size(x,1),1);
  
  Loc = find(abs(x(:,2)) > L);
  sigma(Loc) = sigma_0;

return
