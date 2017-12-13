function L = assemLieAdCochain_1f(Mesh, FHandle,VHandle, JHandle, QuadRule_1D,varargin)
% assemCochain_1f assemble Cochain of one forms
%
%   L = ASSEMCochain_1f(MESH,FHandel) .... 
%
%
%   Copyright 2007-2009 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants
% Assign output arguments
 
 nEdges=size(Mesh.Edges,1);
 L = zeros(nEdges,1);
 nGuass = size(QuadRule_1D.w,1);
 
 options = odeset('Events', @meet,'AbsTol',10^-10);
 for i = 1:nEdges
       
        P1 = Mesh.Coordinates(Mesh.Edges(i,1),:);
        P2 = Mesh.Coordinates(Mesh.Edges(i,2),:);
        % Compute midpoints of all edges
        x = ones(nGuass,1)*P1+QuadRule_1D.x*(P2-P1);
        dS = ones(nGuass,1)*(P2-P1);
        % compute advected solution
        for j=1:nGuass
            if meet(0,x(j,:))
                rhs = @(t,y)[VHandle([y(1) y(2),0])'; ...
                    reshape(reshape(JHandle([y(1),y(2)]),2,2)*reshape(y(3:6),2,2),4,1)];
                y0 = [x(j,:)'; 1 ; 0; 0; 1];
                [t,y,tm,xm]=ode45(rhs,[0,-40],y0,options);
            else
                xm(1:2)=x(j,1:2);
                xm(3:6)= [ 1 0  0 1];
            end
            Fval(j,:) = FHandle(xm(1:2),0,0)*reshape(xm(3:6),2,2);
        end
        L(i) = sum(QuadRule_1D.w.*sum(Fval.*dS,2));
 end
 
    function [fval,isterminal,direction] = meet(t, y)
        fval = 1-max([(y(1)+1)<10^-8, (y(2)+1)<10^-8]);
        isterminal = 1;
        direction = 0;
        return
        


  