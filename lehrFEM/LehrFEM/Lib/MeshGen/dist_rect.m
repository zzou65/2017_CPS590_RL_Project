function dist = dist_rect(x,x0,a,b)
% DIST_RECT Signed distance for a rectangle.
%
%   DIST = DIST_RECT(X,A,B) computes the value of the signed distance
%   function for a rectangle with side lengths A and B and lower left
%   vertex given by X0.
%
%   Example:
%
%   dist = dist_rect(x,x0,a,b);
  
%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  dist = -min(min(min(-x0(2)+x(:,2),x0(2)+b-x(:,2)),-x0(1)+x(:,1)),x0(1)+a-x(:,1));
  
return