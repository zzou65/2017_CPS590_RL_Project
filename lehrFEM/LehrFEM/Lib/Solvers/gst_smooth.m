function x = gst_smooth(x,A,b)
%GST_SMOOTH transposed Gauss-Seidel smoother
%   
%   X = GST_SMOOTH(X,A,B) calculates one Gauss-Seidel iteration for A*X=B.
%   This function can be used as a smoother in MG.
%
%   See also mg, mg_smooth.

%   Copyright 2006-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

U = triu(A);
x = x+U\(b-A*x);