function varargout = psor(A,b,x0,omega,TOL,MAXIT)

% PSOR solves a linear complementary problem(LCP)
%
%   X = PSOR(A,B,X0,OMEGA,TOL,MAXIT) soves the LCP problem given by
%   
%   x'(Ax-b) = 0, x >= 0, Ax-b >= 0 using the projectet SOR algorithm
% 
%   The tolerence TOL is a stopping criterion for the number of psor
%   iteratons. MAXIT specifies the maximum number of iterations
%
%   [X CONV] = PSOR(A,B,OMEGA,TOL,MAXIT) returns a flag that indicates
%   weather the iteration converges or not within the specified number of
%   iterations
%
%   [X CONV ERR] = PSOR(A,B,OMEGA,TOL,MAXIT) returns a vector with the relative 
%   residual ERR at each iteration
%
%   [X CONV ERR IT] = PSOR(A,B,OMEGA,TOL,MAXIT) also returns the number of
%   iteratons that were performed
%
%   Example:
%
%   x = psor(A,b,xo,1.5,1e-5,500)

%   Copyright 2006-2006 Kari Borset
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland  

n=length(b);
x = x0;
% Initialize algorithm
for i = 1:n
    x(i) = max(0,x(i)+omega*(b(i)-A(i,:)*x)/A(i,i));
end

err = 2*TOL;
it = 0;
ERR = [];
while (err>TOL) && (it <= MAXIT),
    it = it+1; x0 = x;
    for i = 1:n
        x(i) = max(0,x(i)+omega*(b(i)-A(i,:)*x)/A(i,i));
    end
    err = norm(x-x0)/norm(b);
    ERR = [ERR err];
end

conv = 1;
if it == MAXIT,
    conv = 0;
end
switch(nargout)
    case 1
        varargout{1} = x;
    case 2
        varargout{1} = x;
        varargout{2} = conv;
    case 3
        varargout{1} = x;
        varargout{2} = conv;
        varargout{3} = ERR;
    case 4
        varargout{1} = x;
        varargout{2} = conv;
        varargout{3} = ERR;
        varargout{4} = it;
end
return