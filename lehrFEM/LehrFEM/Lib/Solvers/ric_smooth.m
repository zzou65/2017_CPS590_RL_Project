function x = ric_smooth(x,A,b,r)
%RIC_SMOOTH Richardson iteration
%   
%   X = RIC_SMOOTH(X,A,B,R) calculates one step of the Richardson iteration
%   for A*X=B.  This function can be used as a smoother in MG.
%   The argument R should be the inverse of the largest eigenvalue of A.
%
%   See also mg, mg_smooth.

%   Copyright 2006-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


x = x+r*(b-A*x);