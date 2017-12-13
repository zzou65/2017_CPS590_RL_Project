function sigma = sigma_x_var(x,ElemFlag,sigma_0,L,varargin)
% SIGMA_Y_VAR Absorption profile.
%
%   SIGMA = SIGMA_Y_VAR(X,ELEMFLAG,SIGMA_0,L) computes the value of the
%   parabolic absorption profile in x-direction.
%
%   Example:
%
%   sigma = sigma_y_var(x,0,sigma_0,L);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  sigma = zeros(size(x,1),1);
  
  Loc = find(abs(x(:,2)) > L);
  sigma(Loc) = sigma_0*(abs(x(Loc,2))-L).^2;
  
return
