function dist = dist_tri(x,a,b,c)
% DIST_TRI Signed distance for a triangle.
%
%   DIST = DIST_TRI(X,A,B,C) computes the value of the signed distance
%   function for a triangle with points a,b,c which are counter-clockwise
%
%   Example:
%
%   dist = dist_tri(x,a,b,c);

%   Copyright 2005-2005 Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

   dist = -min(min(((x(:,1)-a(1))*(a(2)-b(2))-(x(:,2)-a(2))*(a(1)-b(1)))/norm(a-b), ...
        ((x(:,1)-b(1))*(b(2)-c(2))-(x(:,2)-b(2))*(b(1)-c(1)))/norm(b-c)), ...
        ((x(:,1)-c(1))*(c(2)-a(2))-(x(:,2)-c(2))*(c(1)-a(1)))/norm(c-a));

  return