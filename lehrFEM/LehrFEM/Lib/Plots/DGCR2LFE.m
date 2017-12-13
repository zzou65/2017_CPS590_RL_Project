function U = DGCR2LFE(U_DGCR,Mesh,varargin)
% DGCR2LFE Transform a DGCR finite element solution into a LFE solution.
%
%   U_LFE = DGCR2LFE(U_DGCR,MESH) generates a solution vector U_LFE with
%   the values of the DGCR solution estimated at the vertices
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh.
%
%   Example:
%
%   U_LFE = DGCR2LFE(U_DGCR,Mesh);

%   Copyright 2006-2006 Kari Borset
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants


%Initialize constants

nElements = size(Mesh.Elements,1);
nVertices = size(Mesh.Coordinates,1);
U = zeros(nVertices,1); 
N = zeros(nVertices,1);

%For all elements assemble contribution to solution

for i=1:nElements,
    Vertices = Mesh.Elements(i,:);
    alfa = U_DGCR([3*(i-1)+1 3*(i-1)+2 3*(i-1)+3]);
    l=1;
    for j=Vertices(:),
        U(j) = U(j) + sum(alfa(:)) - 2*alfa(l);
        N(j) = N(j) + 1;
        l = l + 1;
    end    
end

U = U ./ N;

return