function [QR_bnd, QR_inn] = split_QuadRule(QuadRule)
% SPLIT_QUADRULE(QUADRULE) splits QUADRULE
%
% Split_Quadrule(QuadRule) splits a QuadRule into boundary
% QR_bnd and inner QR_inn contributions
%
%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants

QR_bnd.x=[]; QR_bnd.w=[];
QR_inn.x=[]; QR_inn.w=[];

Pts=size(QuadRule.w,1);

for i=1:Pts
    if (QuadRule.x(i,1)==0 || QuadRule.x(i,2)==0 || sum(QuadRule.x(i,:))==1)
        QR_bnd.x=[QR_bnd.x; QuadRule.x(i,:)];
        QR_bnd.w=[QR_bnd.w; QuadRule.w(i)];
    else
        QR_inn.x=[QR_inn.x; QuadRule.x(i,:)];
        QR_inn.w=[QR_inn.w; QuadRule.w(i)];
    end
end
        
    
        