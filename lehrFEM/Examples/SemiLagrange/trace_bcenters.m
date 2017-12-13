function pbB = trace_bcenters(Mesh, velocity, Jac, h, varargin)
% trace_bcenter calculates the position of pulled-forward barycenters .
%
%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

nCoordinates = size(Mesh.Coordinates,1);
nEdges =size(Mesh.Edges,1);
nElements=size(Mesh.Elements,1);
nE=size(Mesh.Elements,1);

% barycenters= ...
%     [sum(1/3*reshape(Mesh.Coordinates(Mesh.Elements,1),[],3),2) ...
%     sum(1/3*reshape(Mesh.Coordinates(Mesh.Elements,2),[],3),2)];

%solve ODE for quadrature points
%  rhs=@(t,x) [ sum(velocity([ x(1:nE) x(nE+1:2*nE)]).*[ones(nE,1) zeros(nE,1)],2) ; ...
%      sum(velocity([ x(1:nE) x(nE+1:2*nE)]).*[zeros(nE,1) ones(nE,1)],2) ; ...
%      sum(Jac([ x(1:nE) x(nE+1:2*nE)]).*[ x(2*nE+1:3*nE) zeros(nE,1) x(3*nE+1:4*nE) zeros(nE,1)],2); ...
%      sum(Jac([ x(1:nE) x(nE+1:2*nE)]).*[ zeros(nE,1) x(2*nE+1:3*nE) zeros(nE,1) x(3*nE+1:4*nE)],2); ...
%      sum(Jac([ x(1:nE) x(nE+1:2*nE)]).*[ x(4*nE+1:5*nE) zeros(nE,1) x(5*nE+1:6*nE) zeros(nE,1)],2); ...
%      sum(Jac([ x(1:nE) x(nE+1:2*nE)]).*[ zeros(nE,1) x(4*nE+1:5*nE) zeros(nE,1) x(5*nE+1:6*nE)],2)];
%
%  y_init=[barycenters(:,1);...
%      barycenters(:,2);...
%      ones(nElements,1);...
%      zeros(nElements,1);...
%      zeros(nElements,1);...
%      ones(nElements,1)];
%
%  %[t,y] = ode45(rhs,[0,h],y_init);
%  y = y_init'+h*rhs(0,y_init)';
%
%  data=y(end,:);

% Preallocate memory
pbB = zeros(nElements,7);

% loop over all elements and (back)-trace vertices
for i = 1:nElements

    vid = Mesh.Elements(i,:);

    % Extract nodal directions
    % or solve 3 ODEs \dot(y)=v(y), y(0)=a1,
    %                          \dot(y)=v(y), y(0)=a3,
    %                          \dot(y)=v(y), y(0)=a3;
    % and define v1=y(\tau)-a1
    %                 v2=y(\tau)-a2
    %                 v3=y(\tau)-a3

    %     v1=direction(vid(1),:);
    %     v2=direction(vid(2),:);
    %     v3=direction(vid(3),:);
    %
    %     vb=1/3*(v1+v2+v3);

    a1 = Mesh.Coordinates(vid(1),:);
    a2 = Mesh.Coordinates(vid(2),:);
    a3 = Mesh.Coordinates(vid(3),:);

    b=1/3*(a1+a2+a3);

    %     rhs = @(t,x)[velocity([x(1) x(2)])'; reshape(Jac([x(1) x(2)])*[x(3) x(4); x(5) x(6)],4,1)];
    %rhs = @(t,x)[velocity([x(1) x(2)])'; ...
    %    sum(Jac([x(1) x(2)]).*[x(3) 0 x(4) 0],2); ...
    %   sum(Jac([x(1) x(2)]).*[0 x(3) 0 x(4)],2); ...
    %    sum(Jac([x(1) x(2)]).*[x(5) 0 x(6) 0],2); ...
    %    sum(Jac([x(1) x(2)]).*[0 x(5) 0 x(6)],2)];
    %y0 = [b'; 1 ; 0; 0; 1];
    %[t,y] = ode45(rhs,[0,h],y0);
    %y = y0'+h*rhs(0,y0)';
    %dir = y(end,[1 2])-b;
    dir=h*velocity(b);
    pbB(i,[4:7]) = [1 0 0 1]+h*Jac(b);
    pbB(i,[1 2 3]) = trace_bcenter(i,dir,Mesh);
    %pbB(i,[4:7])=Jac(pbB(i,[1,2]));

    % barycenter
    %     b=barycenters(i,:);
    %     dir = data([i nE+i])-b;
    %     pbB(i,[4:7]) = [data(2*nE+i) data(3*nE+i) data(4*nE+i) data(5*nE+i)];
    %     pbB(i,[1 2 3]) = trace_bcenter(i,dir,Mesh);

end

return