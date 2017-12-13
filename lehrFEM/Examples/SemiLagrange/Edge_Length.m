function el = Edge_Length(Mesh,varargin)
% calculates the length of edges.
%
%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%nCoordinates = size(Mesh.Coordinates,1);
nEdges =size(Mesh.Edges,1);
%nElements=size(Mesh.Elements,1);

% Preallocate memory
el = zeros(nEdges,1);

% loop over all edges
for i = 1:nEdges
 
    % endpoints
    epid = Mesh.Edges(i,:);
    
    el(i)=sqrt(norm(Mesh.Coordinates(epid(1),:)-Mesh.Coordinates(epid(2),:)));
end
  