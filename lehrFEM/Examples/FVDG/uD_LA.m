function u = uD_LA(x,BdFlag,t,nu,varargin)
% UD_LA Dirichlet boundary data.
%
%   U = UD_LA(X,BDFLAG,T,NU) computes the value of the Dirichlet boundary
%   data.
%
%   The smaller the real valued parameter NU the steeper the gradient of
%   the exact solution at X(:,1) = -X(:,2). 
%   
%   Example:
%
%   u = uD_LA([0 0],-1,1,0.01);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  u = -(1-x(:,1).^2).^2.*(1-x(:,2).^2).^2.*atan((x(:,1)+x(:,2)+2*(1-t))/nu);
  
return