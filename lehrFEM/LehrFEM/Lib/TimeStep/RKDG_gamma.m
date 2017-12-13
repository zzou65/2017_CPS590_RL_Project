function a=RKDG_gamma(max);
% alpha parameter for RKDG timestepping
%
%
%
%
%
%
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrumd
%   CH-8092 Zurich, Switzerland

switch max
    case 1
        a=[1];
    case 2
        a=[1/2 1/2];
    case 3
        a=[1/3 1/2 1/6];
    case 4
        a=[3/8 1/3 1/4 1/24];
    case 5
        a=[11/30  3/8 1/6 1/12 1/120];
    case 6
        a=[53/144 11/30 3/16 1/18 1/48 1/720];
    case 7 
        a=[103/280 53/144 11/60 3/48 1/72 1/240 1/5040];
    case 8
        a=[2219/5760 103/280 53/288 11/180 1/64 1/360 1/1440 1/40320]; 
end

return