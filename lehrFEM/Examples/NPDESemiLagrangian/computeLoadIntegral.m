function [L,storeElements] = computeLoadIntegral(U0,QuadRule,Vertices,mesh,tau,Convectionhandle,Dhandle)
% computeLoadIntegral Asseble Load Vector Corresponding to previous 
% time step and lagrangian tracking position.
%
%   L = computeLoadIntegral(U0,QuadRule,Vertices,mesh,tau,Convectionhandle,Dhandle) assembles  
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
storeElements = zeros(1,nPts);
x = QuadRule.x*BK+ones(nPts,1)*bK;

%Preallocate memory to U(x-tau*V)
U = zeros(size(QuadRule.x,1),1);

nElements = size(mesh.Elements,1);

xShift = [0 0];

for i = 1:size(x,1)
%Calculate Shift Position of x to (x-tau*V) by Convectionhandle
    Shift = tau*Convectionhandle(x(i,:));
    
%New position of x as xShift
    xShift(1) = x(i,1) - Shift(1);
    xShift(2) = x(i,2) - Shift(2);
    
%Inner Loop of search for xShift in all the mesh triangles
    for j = 1:nElements
        idx = mesh.Elements(j,:);
        VerticesSearch = mesh.Coordinates(idx,:);
        UNode = U0(idx,:);
        [p,SC] = triangleSearch(VerticesSearch,xShift,UNode);
        if SC == 1
            U(i,1) = p;
            storeElements(1,i) = j;
            break;
        end;
    end;
    if SC == 0
%if search fails then assign boundary condition value to U through Dhandle
        U(i,1) = Dhandle(xShift);
        
    end;
end;

%Finally assign Load with the value of U multiplied weight of quadrature
L = U.*QuadRule.w;



