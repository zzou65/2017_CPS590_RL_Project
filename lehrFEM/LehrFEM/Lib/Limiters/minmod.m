function m=minmod(a, b ,c,M)
% modified (weighted) minmod function
% 
%
%  a,b,c input for minmod
%  M: parameter for high order limiting
%
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

m=0;
if (abs(a)<=M) m=a;
else if (abs(sum(sign([a b c])))==3)
        m=sign(a)*min(abs([a b c]));
    end
end
return;