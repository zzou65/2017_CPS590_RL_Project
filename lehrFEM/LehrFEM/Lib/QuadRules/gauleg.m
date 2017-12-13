function QuadRule = gauleg(a,b,n)
% GL_RULE 1D Gauss-Legendre quadrature rule.
%
%   QUADRULE = GAULEG(a,b,n) computes the N-point Gauss-Legendre
%   quadrature rule on the interval [a,b] up to mashine precision. 
%
%   Note that all quadrature rules obtained from GAULEG are of order 2*N-1.
%
%   The struct QUADRULE contains the following fields:
%    W N-by-1 matrix specifying the weights of the quadrature rule.
%    X N-by-1 matrix specifying the abscissae of the quadrature rule.
%   
%   Example:
%
%   QuadRule = gauleg(0,1,10);

if (n==1), x = 0; w = 2;
else
 c = zeros(n-1,1);
 for i=1:(n-1), c(i)=i/sqrt(4*i*i-1); end
 J=diag(c,-1)+diag(c,1); [ev,ew]=eig(J);
 x=diag(ew); w=(2*(ev(1,:).*ev(1,:)))';
end

% The above rule computed a quadrature rule for the reference interval. 
% The reference interval is always [-1,1], and if [a,b] != [-1,1] then we 
% re-map the abscissas (x) and rescale the weights.

xm = (b+a)/2; % midpoint
xl = (b-a)/2; % area of the requested interval / area of the reference interval

x = xm + xl*x;
w = w * xl;

QuadRule.w = w(:);
QuadRule.x = x(:);
