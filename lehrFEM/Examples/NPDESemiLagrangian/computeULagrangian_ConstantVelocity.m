function [L] = computeULagrangian_ConstantVelocity(U0,Vertices,mesh,tau,Convectionhandle,Dhandle,StoreElements)
% computeULagrangian_ConstantVelocity Assemble Load Vector Corresponding to previous 
% time step and lagrangian tracking position for a stationary convection speed
%
%   L = computeLoadIntegral_ConstantVelocity(U0,QuadRule,Vertices,mesh,tau,Convectionhandle,Dhandle) assembles  
%   Load Vector based on Convectionhandle for IU(x-tau*v)
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%   by Rajdeep Deb
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%Preallocate memory to U(x-tau*V)
U = zeros(3,1);

xShift = [0 0];

for i = 1:3
%Calculate Shift Position of x to (x-tau*V) by Convectionhandle
    Shift = tau*Convectionhandle(Vertices(i,:));
    
%New position of x as xShift
    xShift(1) = Vertices(i,1) - Shift(1);
    xShift(2) = Vertices(i,2) - Shift(2);
    
    ElementLocation = StoreElements(1,i);
    
    if (ElementLocation > 0)
        idx = mesh.Elements(ElementLocation,:);
        Vertices2 = mesh.Coordinates(idx,:);
        UNode = U0(idx,:);
    
        P1 = Vertices2(1,:);
        P2 = Vertices2(2,:);
        P3 = Vertices2(3,:);

        bK = [P2 - P1 ; P3 - P1];

        PT = ((xShift - P1)*inv(bK));
        N = shap_LFE(PT);
        U(i,1) = UNode(1,1)*N(1,1) + UNode(2,1)*N(1,2) + UNode(3,1)*N(1,3);
    else
        U(i,1) = Dhandle(xShift);
    end
end

%Finally assign Load with the value of IU(x-tau*v)
L = U;
