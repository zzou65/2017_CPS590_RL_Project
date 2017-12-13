function x = gs_smooth(x,A,b)
%GS_SMOOTH Gauss-Seidel smoother
%   
%   X = GS_SMOOTH(X,A,B) calculates one Gauss-Seidel iteration for A*X=B.
%   This function can be used as a smoother in MG.
%
%   See also mg, mg_smooth.

%   Copyright 2006-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

L = tril(A);
x = x+L\(b-A*x);