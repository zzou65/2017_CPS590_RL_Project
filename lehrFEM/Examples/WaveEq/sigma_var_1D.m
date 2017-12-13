function sigma = sigma_const_1D(x,sigma_0,L,varargin)
% SIGMA_VAR_1D Absporption profile for 1D PML.
%
%   SIGMA = SIGMA_VAR_1D(X,SIGMA_0,L) computes the function value of the
%   parabolic absorption profile for a right hand side 1D PML layer
%   starting at L. 
%
%   SIGMA_0 is the scaling parameter of the absorption profile.
%
%   Example:
%
%   sigma = sigma_var_1D(x,sigma_0,L);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memroy
  
  sigma = zeros(size(x,1),1);
  
  % Compute function value
  
  loc = find(abs(x) > L);
  sigma(loc) = sigma_0*(abs(x(loc))-L).^2;
  
return