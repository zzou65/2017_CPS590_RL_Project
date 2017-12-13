function f = f_Ball(x,varargin)
% F_BALL Right hand side source.
%
%   F = F_BALL(X) computes the value of the right hand side source term of
%   the Laplace equation for the function 
%
%     u(r,theta) = sin(pi*x1)*sinh(*x2)
%
%   on the unit ball.
%
%   Example:
%
%   f = f_Ball([0 0]);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Compute right hand side source term

 f=zeros(size(x,1),1);
  
return