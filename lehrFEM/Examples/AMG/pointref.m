function eta = pointref(U,Mesh)
%POINTREF refinement towards 0

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  numelem = size(Mesh.Elements,1);

  eta = zeros(numelem,1);

  pnt = find(sum(abs(Mesh.Coordinates),2)==0);

  eta(Mesh.AdjElements(pnt,1:Mesh.nAdjElements(pnt))) = 1;

return
  