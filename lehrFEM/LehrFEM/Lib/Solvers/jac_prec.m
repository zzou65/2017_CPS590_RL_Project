function x = jac_prex(x,A,varargin)
% JAC_PREC Jacobian preconditioner.
%
%   X = JAC_PREC(X,A) computes the value of the Jacobian preconditioner for 
%   the matrix D and the right-hand side x.
%
%   Example:
%
%   x = jac_prec(x,A);

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

    D = diag(A);
    x = x./D;
    
return