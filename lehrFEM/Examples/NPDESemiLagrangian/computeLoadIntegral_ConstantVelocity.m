function [L] = computeLoadIntegral_ConstantVelocity(U0,QuadRule,Vertices,mesh,tau,Convectionhandle,Dhandle,StoreElements)
% computeLoadIntegral_ConstantVelocity Assemble Load Vector Corresponding to previous 
% time step and lagrangian tracking position for a stationary convection speed
%
%   L = computeLoadIntegral_ConstantVelocity(U0,QuadRule,Vertices,mesh,tau,Convectionhandle,Dhandle) assembles  
%   Load Vector based on Convectionhandle
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

%Transformation of x from [0 0;1 0:0 1] Triangle domain into P Triangle
%domain

P1 = Vertices(1,:);
P2 = Vertices(2,:);
P3 = Vertices(3,:);

BK = [ P2 - P1 ; P3 - P1];
bK = P1;
nPts = size(QuadRule.w,1);
x = QuadRule.x*BK+ones(nPts,1)*bK;

%Preallocate memory to U(x-tau*V)
U = zeros(size(QuadRule.x,1),1);

xShift = [0 0];

for i = 1:size(x,1)
%Calculate Shift Position of x to (x-tau*V) by Convectionhandle
    Shift = tau*Convectionhandle(x(i,:));
    
%New position of x as xShift
    xShift(1) = x(i,1) - Shift(1);
    xShift(2) = x(i,2) - Shift(2);
    
    ElementLocation = StoreElements(1,i);
    
    if (ElementLocation > 0)
        idx = mesh.Elements(ElementLocation,:);
        Vertices = mesh.Coordinates(idx,:);
        UNode = U0(idx,:);
    
        P1 = Vertices(1,:);
        P2 = Vertices(2,:);
        P3 = Vertices(3,:);

        bK = [P2 - P1 ; P3 - P1];

        PT = ((xShift - P1)*inv(bK));
        N = shap_LFE(PT);
        U(i,1) = UNode(1,1)*N(1,1) + UNode(2,1)*N(1,2) + UNode(3,1)*N(1,3);
    else
        U(i,1) = Dhandle(xShift);
    end
end

%Finally assign Load with the value of U multiplied weight of quadrature
L = U.*QuadRule.w;
