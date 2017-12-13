function x = mass_prec(x,M)
% MASS_PREC mass preconditioner.
%
%   X = MASS_PREC(X,A) computes the value of the MASS preconditioner for 
%   the matrix M and the right-hand side x.
%
%   Example:
%
%   x = mass_prec(x,M);

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

    x = M\x;
    
return
