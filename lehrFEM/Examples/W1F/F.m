function value = F(x,varargin)
% implementation of the right hand side for Script_W1F_Conv

%   Copyright 2009-2009 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

curlu = pi*(cos(pi*x(:,1))-cos(pi*x(:,2)));

curlcurlu = [pi^2*sin(pi*x(:,2))  pi^2*sin(pi*x(:,1))];
wcrocurlu = [x(:,2).*curlu  -2*x(:,1).*curlu];
u = [sin(pi*x(:,2))  sin(pi*x(:,1))];

value = curlcurlu + wcrocurlu + u;
    