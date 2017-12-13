function pbB = trace_bcenters(Mesh, direction, varargin)
% trace_bcenter calculates the position of pulled-forward barycenters .
%
%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

nCoordinates = size(Mesh.Coordinates,1);
nEdges =size(Mesh.Edges,1);
nElements=size(Mesh.Elements,1);

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

    % barycenter
    b=1/3*(a1+a2+a3);
        
    vb=direction(b);

    pbB(i,[1 2 3])=trace_bcenter(i,vb,Mesh);
    
end
    
return