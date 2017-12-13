function varargout = gradFunc

% GRADFUNC computes the optimal grading function for the minimal surface
% problem with exact solution u = sqrt(cosh(y)^2-x^2) and saves the result
% in gradFunc.mat
%
%   [X G] = GRADFUNC gives the function handles to the uniformly
%   distributet X and the graded G. The grading goes from 0 to 1.
%
%   Example
%
%   [x g] = gradFunc;

% Copyright 2006-2006 Kari Borset
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

clc

TOL = 1e-15;

N = 5000;
h = 1/(N+1);
g = zeros(N+1,1);
g(end) = 1;
g(1)=1;
x = linspace(1,0,N+1)';

cu=-1;
cl = -2;

while abs(g(end))>TOL,

    c = (cu+cl)/2;
    odefun = @(g)c/abs(uxx(g))^(2/3);
    
    for i=1:N,
        g(i+1)=g(i)+h*odefun(g(i));
    end
    if g(end)>0,
        cu = c;
    else
        cl = c;
    end

end

if nargout == 2,
    varargout{1} = x;
    varargout{2} = g;
end

o = [x g];

save gradFunc.mat o

function f=uxx(x)

if abs(x) < 1e-10,
    f = 1e15;
else
    root = sqrt(1-(1-x)^2);
    f = (-root+((1-x)^2)/root)/root;
    
end

return