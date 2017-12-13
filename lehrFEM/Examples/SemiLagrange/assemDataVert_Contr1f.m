function A = assemDataVert_Contr1f(Mesh, vHandle, varargin)
%  assemDataVert_Contr1f calculates data needed in assemMatContr_1FVerUP for assembling 
%  Contraction based on upwind quadrature at vertices

%   Copyright 2007-2009 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants
TOL=eps;

nCoordinates = size(Mesh.Coordinates,1);
nEdges =size(Mesh.Edges,1);
nElements=size(Mesh.Elements,1);

% Preallocate memory
A = zeros(nCoordinates,2);

% look for upwind nodes on each triangle
loc=1:9;
for i = 1:nElements

    % Vertices
    vid = Mesh.Elements(i,:);
    a1 = Mesh.Coordinates(vid(1),:);
    a2 = Mesh.Coordinates(vid(2),:);
    a3 = Mesh.Coordinates(vid(3),:);

    % Compute element mapping

    bK = a1;
    BK = [a2-bK; ...
        a3-bK];
    det_BK = abs(det(BK));
    inv_BK = inv(BK);

    % Extract nodal vectors
    v1=-vHandle(a1);
    v2=-vHandle(a2);
    v3=-vHandle(a3);

    % Compute decision variables Theta(T,

    % for first vertex
    y = a1 + v1;
    yhat = (y-bK)*inv_BK;

    if(yhat(1) >= 0 && yhat(2) >= 0)
        A(vid(1),:)=[i 2/det_BK];
    end

    % for second vertex

    y = a2 + v2;
    yhat = (y-bK)*inv_BK;

    if(yhat(2) >= 0 && sum(yhat) <= 1)
        A(vid(2),:)=[i 2/det_BK];
    end

    % for third index

    y = a3 + v3;
    yhat = (y-bK)*inv_BK;

    if(yhat(1) >= 0 && sum(yhat) <= 1)
        A(vid(3),:)=[i 2/det_BK];
    end
end


return
