function pbV = trace_vertices(Mesh, direction, varargin)
% trace_vertices calculates the position of pulled-forward vertices .
%
%   Mesh
%   direction: array of displacement for each vertex
%   pbV:  [ x-Coord,y_Coord, Elementid]
%
%
%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

nCoordinates = size(Mesh.Coordinates,1);
nEdges =size(Mesh.Edges,1);
nElements=size(Mesh.Elements,1);

% Preallocate memory

pbV = zeros(nCoordinates,3);
%pbV(:,[1,2])=Mesh.Coordinates;

% loop over all elements and (back)-trace vertices
for i = 1:nElements
    
    % Vertices
    vid = Mesh.Elements(i,:);
    a1 = Mesh.Coordinates(vid(1),:);
    a2 = Mesh.Coordinates(vid(2),:);
    a3 = Mesh.Coordinates(vid(3),:);

    % Extract global edge numbers
    eidx = [Mesh.Vert2Edge(vid(2),vid(3)) ...
        Mesh.Vert2Edge(vid(3),vid(1)) ...
        Mesh.Vert2Edge(vid(1),vid(2))];

    % Compute element mapping
    bK = a1;
    BK = [a2-bK; ...
        a3-bK];
    det_BK = abs(det(BK));
    inv_BK = inv(BK);

    % gradient of shape functions
    gradN1=[-1,-1]*inv_BK';
    gradN2=[1,0]*inv_BK';
    gradN3=[0,1]*inv_BK';

    % Extract nodal directions
    % or solve 3 ODEs \dot(y)=v(y), y(0)=a1, 
    %                          \dot(y)=v(y), y(0)=a3,
    %                          \dot(y)=v(y), y(0)=a3;
    % and define v1=y(\tau)-a1
    %                 v2=y(\tau)-a2
    %                 v3=y(\tau)-a3

    v1=direction(vid(1),:);
    v2=direction(vid(2),:);
    v3=direction(vid(3),:);

    % for first vertex
    if norm(v1)-eps>0
        y = a1 + v1;
        yhat = (y-bK)*inv_BK;
        % check if element i intersects with extrusion of first vertex
        if(yhat(1) >= 0 && yhat(2) >= 0)
            % start tracing
            pbV(vid(1),:)=trace_vertex(a1,a2,a3,vid(2),vid(3),gradN1,i,v1,Mesh);
        end
    else
        pbV(vid(1),:)=[a1,i];
    end

    % for second vertex
    if norm(v2)-eps>0
        y = a2 + v2;
        yhat = (y-bK)*inv_BK;
        % check if element i intersects with extrusion of second vertex
        if(yhat(2) >= 0 && sum(yhat) <= 1)
            %start tracing
            pbV(vid(2),:)=trace_vertex(a2,a1,a3,vid(1),vid(3),gradN2,i,v2,Mesh);
            %trace_vertex(vertex,vertex_r,vertex_q,r_id,q_id,gradN,Element,velo,Mesh)
        end
    else
        pbV(vid(2),:)=[a2,i];
    end

    % for third index
    if norm(v3)>eps
        y = a3 + v3;
        yhat = (y-bK)*inv_BK;
        % check if element i intersects with extrusion of second vertex
        if(yhat(1) >= 0 && sum(yhat) <= 1)
            % start traceing
            pbV(vid(3),:)=trace_vertex(a3,a1,a2,vid(1),vid(2),gradN3,i,v3,Mesh);
            %trace_vertex(vertex,vertex_r,vertex_q,r_id,q_id,gradN,Element,
            %velo,Mesh)
        end
    else
        pbV(vid(3),:)=[a3,i];
    end
end

% Assign output arguments
return