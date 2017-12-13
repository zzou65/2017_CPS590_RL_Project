function c=innerproduct(F1Handle,F2Handle, xmin,xmax,ymin,ymax);
% innerproduct calculates the innerproduct on square xmin,xmax,ymin,ymax
%
%   Example:
%   F1Handle=@(x) x(:,1)+x(:,2)
%   F2Handle=@(x) x(:,1);
%
%   c=innerproduct(F1Handle,F2Handle,-1,1,-1,1)
%
%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

Product=@(x)sum(F1Handle(x).*F2Handle(x),2);
integrand=@(x,y)Product([(size(x,2)==1)*ones(size(y,2),1)*x'+(size(x,2)~=1)*x',...
                                     (size(y,2)==1)*ones(size(x,2),1)*y'+(size(y,2)~=1)*y']);
c=dblquad(integrand,xmin,xmax,ymin,ymax);
end