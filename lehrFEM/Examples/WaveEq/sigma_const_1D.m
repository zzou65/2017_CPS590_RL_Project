function sigma = sigma_const_1D(x,sigma_0,L,varargin)
% SIGMA_CONST_1D Absporption profile for 1D PML.
%
%   SIGMA = SIGMA_CONST_1D(X,SIGMA_0,L) computes the function value of the
%   constant absorption profile for a 1D PML layer. 
%
%   SIGMA_0 is the scaling parameter of the absorption profile.
%   
%   Example:
%
%   sigma = sigma_const_1D(x,sigma_0,L);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Preallocate memory
  
  sigma = zeros(size(x,1),1);
  
  % Compute function value
  
  loc = find(abs(x) > L);
  sigma(loc) = sigma_0;
  
return