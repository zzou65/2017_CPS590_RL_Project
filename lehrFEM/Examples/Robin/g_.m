function g_N = g_(x,BdFlag,varargin)
% G_ parameter for robin boundary condition.
%
%   G_N = sig_(X) computes the value of sigma for the Robin boundary data for
%   the function 
%
%
%   on the domain according to the boundary flag BDFLAG.
%
%   Example:
%
%   g_N = g_([0 0]);

%   Copyright 2005-2006 Patrick Meury & Kah Ling Sia & HOlger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


  % Compute sigma data
  
   % i.)
   %g_N = ones(size(x,1),1);  
   % ii)
   g_N= -1/2.*ones(size(x,1),1);
   
return