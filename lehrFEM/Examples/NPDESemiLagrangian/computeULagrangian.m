function [L,storeElements] = computeULagrangian(U0,Vertices,mesh,tau,Convectionhandle,Dhandle)
% computeULagrangian Assemble Load Vector Corresponding to previous 
% time step and lagrangian tracking position for IU(x-v*tau).
%
%   [L,storeElements] = computeLoadIntegral(U0,QuadRule,Vertices,mesh,tau,Convectionhandle,Dhandle) assembles  
%   Load Vector based on Convectionhandle and storeElements is obtained to
%   be used for time independent velocity field
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

storeElements = zeros(1,3);

%Preallocate memory to U(x-tau*V)
U = zeros(3,1);

nElements = size(mesh.Elements,1);

xShift = [0 0];

for i = 1:3
%Calculate Shift Position of x to (x-tau*V) by Convectionhandle
    Shift = tau*Convectionhandle(Vertices(i,:));
    
%New position of x as xShift
    xShift(1) = Vertices(i,1) - Shift(1);
    xShift(2) = Vertices(i,2) - Shift(2);
    
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

%Finally assign Load with the value of IU(x-tau*v)
L = U;
