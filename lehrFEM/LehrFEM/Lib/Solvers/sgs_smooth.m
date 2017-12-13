function x = sgs_smooth(x,A,b)
%SGS_SMOOTH symmetric Gauss-Seidel smoother
%   
%   X = SGS_SMOOTH(X,A,B) calculates one symmetric Gauss-Seidel iteration
%   for A*X=B.  This function can be used as a smoother in MG.
%
%   See also mg, mg_smooth.

%   Copyright 2006-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

L = tril(A);
U = triu(A);
x = x+L\(b-A*x);
x = x+U\(b-A*x);