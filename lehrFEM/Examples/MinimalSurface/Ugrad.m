function f = Ugrad(x,varargin)

% UGRAD computes the gradient of U = sqrt(cosh(y).^2 - x.^2)
%
%   X is a Nx2 vector with x and y coordinates of the N points 
%
%   Example
%
%   f = Ugrad([.5 .5 ; .3 .3]);

% Copyright 2006-2006 Kari Borset
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

n = size(x,1);
f = zeros(n,2);

for i = 1:n,
    coord = x(i,:);
    Ue = sqrt(cosh(coord(2)).^2 - coord(1).^2);
    f(i,:) = 1/Ue * [-coord(1) .5*sinh(2*coord(2))];
end