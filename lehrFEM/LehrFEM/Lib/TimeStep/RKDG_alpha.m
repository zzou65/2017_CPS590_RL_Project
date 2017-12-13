function a=RKDG_alpha(max);
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
        a=[1 0; 1/2 1/2];
    case 3
        a=[1 0 0; 3/4 1/4 0; 1/3 0 2/3];
end

return