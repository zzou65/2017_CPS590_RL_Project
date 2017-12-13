function [Loc1, Loc2] = connect_period(Mesh, BdFlag1, BdFlag2, varargin)
%   find corresponding DOFS for periodic boundary conditions
%
%
%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants

% Extract edges
  
Loc = get_BdEdges(Mesh);
Loc1 = Loc(Mesh.BdFlags(Loc) == BdFlag1);
Loc2 = Loc(Mesh.BdFlags(Loc) == BdFlag2);

ID1=Mesh.Edges(Loc1,:);
ID2=Mesh.Edges(Loc2,:);

MidPoints1_y=1/2*(Mesh.Coordinates(ID1(:,1), 2 ) + Mesh.Coordinates(ID1(:,2), 2 ));
MidPoints2_y=1/2*(Mesh.Coordinates(ID2(:,1), 2 ) + Mesh.Coordinates(ID2(:,2), 2 ));

[MidPoints1_y, id1] = sort(MidPoints1_y);
[MidPoints2_y, id2] = sort(MidPoints2_y);

Loc1 = Loc1(id1);
Loc2 = Loc2(id2);

return