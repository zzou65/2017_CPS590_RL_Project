function j = divj_LShap_Hdiv(x,varargin)
% divj_LShap_Hdiv Exact solution for H div semi-norm.
%
%   UEX = divj_LShap_Hdiv(X)  the laplace of the function
%
%      u(r,theta) = r^2/3*sin(2/3*theta)
%
%   at the point X.
%
%   Example:
%
%   uex = divj_LShap_Hdiv([0 0])

%   Copyright 2005-2006 Patrick Meury & Kah Ling Sia & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 
  
  j=zeros(size(x,1),1);
  
return