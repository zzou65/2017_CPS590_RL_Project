function dist = dist_circ(x,c,r)
% DIST_CIRC Signed distance for a circle.
%
%   DIST = DIST_CIRC(X,C,R) computes the value of the signed distance
%   function of the circle with radius R and center C.
%
%   Example:
%
%   dist = dist_circ(x,c,r);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  dist = sqrt((x(:,1)-c(1)).^2+(x(:,2)-c(2)).^2)-r;
  
return