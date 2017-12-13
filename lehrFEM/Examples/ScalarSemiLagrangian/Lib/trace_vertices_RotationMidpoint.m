function pbV = trace_vertices_RotationMidpoint(Mesh, h,nsteps_loc, rule, varargin)
%  trace_vertices calculates the position of pulled-forward vertices .
%
%  Copyright 2008-2008 Holger Heumann
%  SAM - Seminar for Applied Mathematics
%  ETH-Zentrum
%  CH-8092 Zurich, Switzerland

nCoordinates = size(Mesh.Coordinates,1);
nEdges =size(Mesh.Edges,1);
nElements=size(Mesh.Elements,1);

% calcute trajectories via local timestepping

nC=nCoordinates;
Cx=Mesh.Coordinates(:,1);
Cy=Mesh.Coordinates(:,2);

Cx0=Mesh.Coordinates(:,1);
Cy0=Mesh.Coordinates(:,2);

tau = h/nsteps_loc;

for j=1:nsteps_loc

    % implicit midpoint rule
    if  strcmp(rule,'Midpoint Rule')
        directionX=1./(1+tau^2/4).*((1-tau^2/4).*Cx0+tau*Cy0)-Cx0;
        directionY=-tau/2*(Cx0+1/(1+tau^2/4).*((1-tau^2/4).*Cx0+tau*Cy0));
    end

    % explicit Euler rule
    if strcmp(rule,'Explicit Euler')
        directionX=tau*Cy0;
        directionY=-tau*Cx0;
    end

    % implicit Euler rule
    if strcmp(rule, 'Implicit Euler')
        directionX=1./(1+tau^2).*(Cx0+tau*Cy0)-Cx0;
        directionY=-tau*1./(1+tau^2).*(Cx0+tau*Cy0);
    end

    % Trapez rule
    if strcmp(rule, 'Trapez Rule')
        k1X = Cy0;
        k1Y = - Cx0;
        k2X = Cy0 + tau * k1Y;
        k2Y = - Cx0 - tau * k1X;
        directionX = tau/2*(k1X+k2X);
        directionY = tau/2*(k1Y+k2Y);
    end

    % Explicit Midpoint Rule
    if strcmp(rule, 'Explicit Midpoint Rule')
        k1X = Cy0;
        k1Y = - Cx0;
        k2X = Cy0 + tau/2 * k1Y;
        k2Y = - Cx0 - tau/2 * k1X;
        directionX = tau*(k2X);
        directionY = tau*(k2Y);
    end

    if strcmp(rule,'exact')
        directionX = (cos(-tau)*Cx0-sin(-tau)*Cy0)-Cx0;
        directionY = (sin(-tau)*Cx0+cos(-tau)*Cy0)-Cy0;
    end

    Cx0=Cx0+directionX;
    Cy0=Cy0+directionY;
end

direction=[Cx0-Cx Cy0-Cy];

% Preallocate memory
pbV = [Mesh.Coordinates,zeros(nC,1)];

% loop over all boundary edges and (back)-trace vertices
bndEdges=get_BdEdges(Mesh);
bndVertices= unique(Mesh.Edges(bndEdges,:));

for i =1:size(bndVertices)
    pbV(bndVertices(i),[1 2 3])=[Mesh.Coordinates(bndVertices(i),:)+direction(bndVertices(i),:),0];
end


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
            pbV(vid(1),[1 2 3])=trace_vertex(a1,a2,a3,vid(2),vid(3),gradN1,i,v1,Mesh);
        end
    else
        pbV(vid(1),[1 2 3])=[a1,i];
    end

    % for second vertex
    if norm(v2)-eps>0
        y = a2 + v2;
        yhat = (y-bK)*inv_BK;
        % check if element i intersects with extrusion of second vertex
        if(yhat(2) >= 0 && sum(yhat) <= 1)
            %start tracing
            pbV(vid(2),[1 2 3])=trace_vertex(a2,a1,a3,vid(1),vid(3),gradN2,i,v2,Mesh);
            %trace_vertex(vertex,vertex_r,vertex_q,r_id,q_id,gradN,Element,velo,Mesh)
        end
    else
        pbV(vid(2),[1 2 3])=[a2,i];
    end

    % for third index
    if norm(v3)>eps
        y = a3 + v3;
        yhat = (y-bK)*inv_BK;
        % check if element i intersects with extrusion of second vertex
        if(yhat(1) >= 0 && sum(yhat) <= 1)
            % start traceing
            pbV(vid(3),[1 2 3])=trace_vertex(a3,a1,a2,vid(1),vid(2),gradN3,i,v3,Mesh);
            %trace_vertex(vertex,vertex_r,vertex_q,r_id,q_id,gradN,Element,
            %velo,Mesh)
        end
    else
        pbV(vid(3),[1 2 3])=[a3,i];
    end
end

% Assign output arguments
return
