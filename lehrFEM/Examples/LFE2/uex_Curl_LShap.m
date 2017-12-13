function uex = uex_Curl_LShap(x,varargin)
% UEX_CURL_LSHAP Exact solution for H1 semi-norm.
%
%   UEX = UEX_CURL_LSHAP(X) computes curl of the gradient of the function
%
%      u(r,theta) = r^2/3*sin(2/3*theta)
%
%   at the point X.
%
%   Example:
%
%   uex = uex_Curl_LShap([0 0]);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  % Compute the function value

  uex = zeros(size(x,1),1);
  
return