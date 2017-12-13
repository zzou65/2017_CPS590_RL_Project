function QuadRule=gaulob(a,b,N,tol)
% GL_RULE 1D Legendre-Gauss-Lobatto quadrature rule.
%
%   QUADRULE = GAULOB(A,B,N,TOL) computes the N-point Legendre-Gauss
%   Lobatto quadrature rule on the interval [A,B] up to the prescribed tolerance
%   TOL. If no tolerance is prescribed GAULEG uses the machine precision
%   EPS.
%
%   Note that all quadrature rules obtained from GAULOB are of order 2*N-3.
%
%   The struct QUADRULE contains the following fields:
%    W N-by-1 matrix specifying the weights of the quadrature rule.
%    X N-by-1 matrix specifying the abscissae of the quadrature rule.
%   
%   Example:
%
%   QuadRule = gaulob(0,1,10,1e-6);

%   Copyright 2005-2008 Patrick Meury, Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Check for right number of arguments 

if(nargin == 3)
    tol = eps;          
end

xm=(a+b)/2;
xl=(b-a)/2;

% Truncation + 1
N1=N+1;

% Use the Chebyshev-Gauss-Lobatto nodes as the first guess
x=cos(pi*(0:N-1)/(N-1))';

% The Legendre Vandermonde Matrix
P=zeros(N,N);

% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and 
% update x using the Newton-Raphson method.

xold=2;

while max(abs(x-xold))>tol

    xold=x;
	        
    P(:,1)=1;    P(:,2)=x;
    
    for k=2:N
	P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
		    
    x=xold-( x.*P(:,N)-P(:,N-1) )./( N1*P(:,N));
				                 
end
QuadRule.x=xm+xl*x;
QuadRule.w=2*xl./((N-1)*N*P(:,N).^2);
