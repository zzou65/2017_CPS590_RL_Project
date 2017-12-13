function v= velo(x,varargin)
% F_LSHAP Right hand side source.
%
%   F = F_LSHAP(X) computes the value of the right hand side source term of
%   the Laplace equation for the function 
%
%     u(r,theta) = r^2/3*sin(2/3*theta)
%
%   on the L-shaped domain.
%
%   Example:
%
%   f = f_LShap([0 0]);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Compute right hand side source term
  
  v = 1/sqrt(10)*[1,3];
  %v = 1/sqrt(1)*[1,0];
return